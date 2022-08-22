#include "node.hpp"

void Node::add_child(Node *child) {
    children_.push_back(child);
    child->parent_ = this;
}
size_t Node::get_id() {
    return id_;
}
const std::vector<Node *> &Node::get_children() {
    return children_;
}
