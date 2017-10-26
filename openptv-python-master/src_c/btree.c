/*  Implementation of btree.h - a binary tree holding a general object at each
    node.
    
    References:
    [1] http://en.wikipedia.org/wiki/Binary_tree
    [2] http://en.wikipedia.org/wiki/Tree_traversal#Depth-first
*/

#include <stdlib.h>
#include "btree.h"

/*  Insert an item to the tree using the item's sort value. The item becomes
    owned by the tree and is freed when the tree is destroyed.
    
    Arguments:
    btree_node **tree - points to the pointer holding the tree root. That's
        because if that pointer is NULL we allocate a new btree_node.
    void *val - a pointe rto the value stored in the node. The tree takes
        ownership of the memory handle.
    double sort_value - when the tree is traversed in-order [2], the items
        will be returned in increasing sort-value order.
*/
void btree_insert(btree_node **tree, void *val, double sort_value) {
    if( *tree == NULL ) {
        *tree = (btree_node *) malloc( sizeof(btree_node) );
        (*tree)->sort_value = sort_value;
        (*tree)->val = val;
        
        /* initialize the children to null */
        (*tree)->left = NULL;
        (*tree)->right = NULL;
    } else if(sort_value <= (*tree)->sort_value) {
        btree_insert(&((*tree)->left), val, sort_value);
    } else if(sort_value > (*tree)->sort_value) {
        btree_insert(&((*tree)->right), val, sort_value);
    }
}

/*  btree_free() recursively destroys a tree and frees its memory and the memory
    held by node items.
    
    Arguments:
    btree_node *tree - the root element of the tree.
*/
void btree_free(btree_node *tree) {
    if (tree == NULL) return;
    
    btree_free(tree->left);
    btree_free(tree->right);
    free(tree->val);
    free(tree);
}

