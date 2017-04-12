#ifndef BIG_UNORDERED_MAP_H_
#define BIG_UNORDERED_MAP_H_

#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#include <cstddef>
#include <functional>
#include <limits>
#include <list>
#include <unordered_map>
#include <vector>

// An unordered map for storing large amount of data over several MPI nodes.
// K: key, V: value, H: hasher.
// Maximum size: 2^64 - 1.
// Maximum async calls between complete_async: 2^64 - 1.
template <class K, class V, class H>
class BigUnorderedMap {
 public:
  BigUnorderedMap(
      boost::mpi::communicator& world,
      const std::pair<K, V>& skeleton,
      const int total_buf_size = 2000);
  void reserve(const unsigned long long);
  unsigned long long bucket_count(const int root = 0) const;
  unsigned long long size(const int root = 0) const;
  void async_inc(const K&, const V&);
  void complete_async_incs();  // Must be called to complete async inc reqs.
  const std::unordered_map<K, V, H>& get_local_map() { return local_map; }
  static void all_combine(
      std::unordered_map<K, V, H>&, boost::mpi::communicator&);

 private:
  void complete_async_inc_trunk();  // Used by async inc when buf full.
  void set_skeleton(const std::pair<K, V>&);  // Called at the beginning.
  void reset_cnts();  // Reset send and recv counts.
  void set_node_buckets();
  int proc_id;
  int n_procs;
  std::unordered_map<K, V, H> local_map;
  int buf_size;
  int send_buf_cnt;
  // buf_send(recv)_content stores content of buf_send(recv)
  // MPI sends content only without structure.
  std::vector<std::pair<K, V>> buf_send;
  std::vector<boost::mpi::content> buf_send_content;
  std::pair<K, V> buf_recv;
  boost::mpi::content buf_recv_content;
  H hasher;
  boost::mpi::communicator* world;
  std::list<boost::mpi::request> reqs;
  std::vector<unsigned long long> send_cnts;
  std::vector<unsigned long long> recv_cnts;
  std::vector<unsigned long long> recv_totals;
  std::vector<unsigned long long> recv_trunk_totals;
  enum { TAG_NODES, TAG_KV, TAG_TRUNK_FINISH, TAG_FINISH };
  int total_node_buckets;
  int local_node_buckets;
  std::vector<int> node_map;
};

template <class K, class V, class H>
BigUnorderedMap<K, V, H>::BigUnorderedMap(
    boost::mpi::communicator& world,
    const std::pair<K, V>& skeleton,
    const int total_buf_size) {
  this->world = &world;
  proc_id = world.rank();
  n_procs = world.size();
  hasher = H();
  buf_size = std::max(total_buf_size / n_procs, 100);
  buf_send.resize(buf_size);
  buf_send_content.resize(buf_size);
  reset_cnts();
  this->set_skeleton(skeleton);
  this->set_node_buckets();
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::reserve(const unsigned long long n_buckets) {
  std::size_t local_buckets = static_cast<std::size_t>(
      n_buckets * 1.0 * local_node_buckets / total_node_buckets + 1);
  printf(
      "Proc %d about to reserve %'lu hash buckets.\n", proc_id, local_buckets);
  local_map.reserve(local_buckets);
  printf(
      "Proc %d reserved %'lu local hash buckets.\n",
      proc_id,
      local_map.bucket_count());
  world->barrier();
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::set_skeleton(const std::pair<K, V>& skeleton) {
  for (int i = 0; i < buf_size; i++) {
    buf_send[i] = skeleton;
    buf_send_content[i] = boost::mpi::get_content(buf_send[i]);
  }
  buf_recv = skeleton;
  buf_recv_content = boost::mpi::get_content(buf_recv);
}

template <class K, class V, class H>
unsigned long long BigUnorderedMap<K, V, H>::bucket_count(
    const int root) const {
  // Only root process gets result.
  unsigned long long local_bucket_count =
      static_cast<unsigned long long>(local_map.bucket_count());
  unsigned long long total_bucket_count = 0;
  reduce(
      *world,
      local_bucket_count,
      total_bucket_count,
      std::plus<unsigned long long>(),
      root);
  return total_bucket_count;
}

template <class K, class V, class H>
unsigned long long BigUnorderedMap<K, V, H>::size(const int root) const {
  // Only root process gets result.
  unsigned long long local_size =
      static_cast<unsigned long long>(local_map.size());
  unsigned long long total_size = 0;
  reduce(*world, local_size, total_size, std::plus<unsigned long long>(), root);
  return total_size;
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::async_inc(const K& key, const V& value) {
  const int target = node_map[hasher(key) % total_node_buckets];

  // Process locally.
  if (target == proc_id) {
    local_map[key] += value;
    send_cnts[proc_id]++;
    return;
  }

  // Send to target asynchronously.
  buf_send[send_buf_cnt].first = key;
  buf_send[send_buf_cnt].second = value;
  reqs.push_front(world->isend(target, TAG_KV, buf_send_content[send_buf_cnt]));
  send_cnts[target]++;
  if (send_cnts[target] == std::numeric_limits<unsigned long long>::max()) {
    printf("Maximum sends reached. Call complete async first.\n");
    exit(1);
  }
  send_buf_cnt++;

  if (send_buf_cnt == buf_size) {
    // When buffer is full, wait and process a trunk.
    // A trunk contains 'buf_size' send-requests from each process.
    for (int i = 0; i < n_procs; i++) {
      if (i == proc_id) continue;
      reqs.push_front(world->isend(i, TAG_TRUNK_FINISH, send_cnts[i]));
    }
    complete_async_inc_trunk();
  }
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::complete_async_inc_trunk() {
  int n_active_procs = 0;
  int trunk_finish_cnt = 0;
  for (int i = 0; i < n_procs; i++) {
    if (recv_cnts[i] < recv_totals[i] && i != proc_id) n_active_procs++;
  }
  while (trunk_finish_cnt < n_active_procs) {
    const auto& status = world->probe();
    const int source = status.source();
    const int tag = status.tag();
    switch (tag) {
      case TAG_KV: {
        world->recv(source, tag, buf_recv_content);
        const K& key = buf_recv.first;
        const V& value = buf_recv.second;
        local_map[key] += value;
        recv_cnts[source]++;
        if (recv_cnts[source] == recv_trunk_totals[source]) {
          trunk_finish_cnt++;
        }
        break;
      }
      case TAG_TRUNK_FINISH: {
        world->recv(source, TAG_TRUNK_FINISH, recv_trunk_totals[source]);
        if (recv_trunk_totals[source] == recv_cnts[source]) {
          trunk_finish_cnt++;
        }
        break;
      }
      case TAG_FINISH: {
        world->recv(source, TAG_FINISH, recv_totals[source]);
        if (recv_totals[source] == recv_cnts[source]) trunk_finish_cnt++;
        break;
      }
    }
  }
  wait_all(reqs.begin(), reqs.end());
  reqs.clear();
  send_buf_cnt = 0;
  recv_trunk_totals.assign(
      n_procs, std::numeric_limits<unsigned long long>::max());
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::complete_async_incs() {
  int n_active_procs = 0;
  for (int i = 0; i < n_procs; i++) {
    if (i == proc_id) continue;
    reqs.push_front(world->isend(i, TAG_FINISH, send_cnts[i]));
    if (recv_cnts[i] < recv_totals[i]) n_active_procs++;
  }
  while (n_active_procs > 0) {
    const auto& status = world->probe();
    const int source = status.source();
    const int tag = status.tag();

    switch (tag) {
      case TAG_KV: {
        world->recv(source, tag, buf_recv_content);
        const K& key = buf_recv.first;
        const V& value = buf_recv.second;
        local_map[key] += value;
        recv_cnts[source]++;
        if (recv_cnts[source] == recv_totals[source]) n_active_procs--;
        break;
      }
      case TAG_TRUNK_FINISH: {
        world->recv(source, TAG_TRUNK_FINISH, recv_trunk_totals[source]);
        break;
      }
      case TAG_FINISH: {
        world->recv(source, TAG_FINISH, recv_totals[source]);
        if (recv_totals[source] == recv_cnts[source]) n_active_procs--;
        break;
      }
    }
  }
  wait_all(reqs.begin(), reqs.end());
  reqs.clear();
  reset_cnts();
  world->barrier();
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::reset_cnts() {
  send_buf_cnt = 0;
  send_cnts.assign(n_procs, 0);
  recv_cnts.assign(n_procs, 0);
  recv_totals.assign(n_procs, std::numeric_limits<unsigned long long>::max());
  recv_trunk_totals.assign(
      n_procs, std::numeric_limits<unsigned long long>::max());
}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::all_combine(
    std::unordered_map<K, V, H>& local_map, boost::mpi::communicator& world) {}

template <class K, class V, class H>
void BigUnorderedMap<K, V, H>::set_node_buckets() {
  // Obtain local node buckets. -1 if not available.
  std::ifstream nodes_file("NODES");
  std::string local_proc_name = boost::mpi::environment::processor_name();
  std::string proc_name_in;
  int proc_mem_in;
  local_node_buckets = 0;
  while (nodes_file >> proc_name_in) {
    nodes_file >> proc_mem_in;
    if (proc_name_in == local_proc_name) {
      local_node_buckets = proc_mem_in;
      break;
    }
  }
  std::list<boost::mpi::request> reqs;
  for (int i = 0; i < n_procs; i++) {
    if (i == proc_id) continue;
    reqs.push_front(world->isend(i, TAG_NODES, local_node_buckets));
  }

  // Receive node buckets info from remote nodes.
  int recv_cnt = 0;
  std::vector<int> node_buckets(n_procs, 0);
  node_buckets[proc_id] = local_node_buckets;
  while (recv_cnt < n_procs - 1) {
    const auto& status = world->probe(boost::mpi::any_source, TAG_NODES);
    recv_cnt++;
    const int source = status.source();
    world->recv(source, TAG_NODES, node_buckets[source]);
  }
  wait_all(reqs.begin(), reqs.end());

  // Setup node map.
  const bool all_unknown =
      std::all_of(node_buckets.begin(), node_buckets.end(), [](const int i) {
        return i == 0;
      });
  if (all_unknown) {
    // No node infomation available, distribute evenly.
    total_node_buckets = n_procs;
    local_node_buckets = 1;
    node_map.reserve(n_procs);
    for (int i = 0; i < n_procs; i++) node_map.push_back(i);
  } else {
    int min_node_buckets = std::numeric_limits<int>::max();
    for (int i = 0; i < n_procs; i++) {
      if (node_buckets[i] > 0 && node_buckets[i] < min_node_buckets) {
        min_node_buckets = node_buckets[i];
      }
    }
    total_node_buckets = 0;
    for (int i = 0; i < n_procs; i++) {
      if (node_buckets[i] == 0) node_buckets[i] = min_node_buckets;
      total_node_buckets += node_buckets[i];
    }
    local_node_buckets = node_buckets[proc_id];
    node_map.reserve(total_node_buckets);
    for (int i = 0; i < n_procs; i++) {
      node_map.insert(node_map.end(), node_buckets[i], i);
    }
  }

  if (proc_id == 0) {
    printf("node_map: ");
    for (const auto i : node_map) printf("%d ", i);
    printf("\n");
  }
  world->barrier();
  printf(
      "Node #%d (%s) stores: %.2f %%\n",
      proc_id,
      local_proc_name.c_str(),
      local_node_buckets * 100.0 / total_node_buckets);
  world->barrier();
  // exit(0);
}

#endif  // BIG_UNORDERED_MAP_H_