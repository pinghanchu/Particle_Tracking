/* Binary tree, useful for sorting etc. */

typedef struct btree_node {
    void *val;
    double sort_value;
    
    struct btree_node *left, *right;
} btree_node;

void btree_insert(btree_node **tree, void *val, double sort_value);
void btree_free(btree_node *tree);
