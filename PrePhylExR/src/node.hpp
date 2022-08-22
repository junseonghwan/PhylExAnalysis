#include <string>
#include <vector>

class Node {
  size_t id_;
  std::vector<Node *> children_;
  Node *parent_ = 0;
public:
  Node(size_t id) : id_(id) {}
  void add_child(Node *child);
  size_t get_id();
  const std::vector<Node *> &get_children();
};

