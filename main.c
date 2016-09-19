/* 
   Meta-LPAm+ community detection algorithm.
   Copyright (C) 2016  Ba-Dung Le <dungleba@gmail.com>   
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

*/

#include <stdio.h>
#include <stdlib.h>
#include <igraph.h>

#define EPSILON 1.e-10

// Meta-LPAm+ reuses some of the data structures and functions that have been implemented for the multi-level modularity optimization algorithm in igraph        

/* Structure storing a community */
typedef struct {
  igraph_integer_t size;           /* Size of the community */
  igraph_real_t weight_inside;     /* Sum of edge weights inside community */
  igraph_real_t weight_all;        /* Sum of edge weights starting/ending
                                      in the community */  
} igraph_i_multilevel_community;

/* Global community list structure */
typedef struct {
  long int communities_no, vertices_no;  /* Number of communities, number of vertices */
  igraph_real_t weight_sum;              /* Sum of edges weight in the whole graph */
  igraph_i_multilevel_community *item;   /* List of communities */
  igraph_vector_t *membership;           /* Community IDs */
  igraph_vector_t *weights;        /* Graph edge weights */
} igraph_i_multilevel_community_list;

/* Computes the modularity of a community partitioning */
igraph_real_t igraph_i_multilevel_community_modularity(const igraph_i_multilevel_community_list *communities) {
  igraph_real_t result = 0;
  long int i;
  igraph_real_t m = communities->weight_sum;
  
  for (i = 0; i < communities->vertices_no; i++) {
    if (communities->item[i].size > 0) {
      result += (communities->item[i].weight_inside - communities->item[i].weight_all*communities->item[i].weight_all/m)/m;
    } 
  }

  return result;
}

// computes the normalized gain in modularity (not exact gain)
igraph_real_t igraph_i_multilevel_community_modularity_gain(const igraph_i_multilevel_community_list *communities, igraph_integer_t community, igraph_integer_t vertex, igraph_real_t weight_all, igraph_real_t weight_inside) {
  return weight_inside - (weight_all * communities->item[(long int)community].weight_all)/communities->weight_sum;
}

// end of reused data structure and functions

struct bin_tree {
    long int to;
    long int edges;
    struct bin_tree * right, * left;
};
typedef struct bin_tree node_edges_to;

struct bin_tree_delta_q {
    double delta_q;
    long int from;    
    long int to; 
    long int edges;    
    struct bin_tree_delta_q * right, * left;
};
typedef struct bin_tree_delta_q node_delta_q;

struct queue{
    double delta_q;
    long int from;
    long int to;    
    long int edges;
};
typedef struct queue node_delta_q_list;

long int max_n_level=1000, n_level=0;
    
void add_edges(node_edges_to **tree, long int to_community, long int val){    
    if(*tree==NULL)
    {
        node_edges_to *temp = NULL;
        temp = (node_edges_to *)malloc(sizeof(node_edges_to));
        temp->left = temp->right = NULL;
        temp->to = to_community;
        temp->edges = val;
        *tree = temp;
        return;        
    }
    
    if(to_community == (*tree)->to)
    {
        (*tree)->edges = (*tree)->edges + val;
    }
    else if(to_community < (*tree)->to)
    {
        add_edges(&(*tree)->left, to_community, val);
    }
    else if(to_community > (*tree)->to)
    {
        add_edges(&(*tree)->right, to_community, val);
    }
}

node_delta_q *find_max(node_delta_q *delta_q){
    if (delta_q==NULL){       
        return NULL;
    }
    if (delta_q->right==NULL){       
        return delta_q;
    }
    return find_max(delta_q->right);    
}

void free_delta_q(node_delta_q **delta_q){
     if (*delta_q==NULL){       
        return;
    }     
     free_delta_q(&(*delta_q)->left);
     free_delta_q(&(*delta_q)->right);
     free(*delta_q);
     IGRAPH_FINALLY_CLEAN(1);
     return;
}

void free_edges_to(node_edges_to **edges_to){
     if (*edges_to==NULL){       
        return;
    }     
     free_edges_to(&(*edges_to)->left);
     free_edges_to(&(*edges_to)->right);
     free(*edges_to);
     IGRAPH_FINALLY_CLEAN(1);
     return;
}

void add_delta_q(node_delta_q **delta_q, double val, long int cfrom, long int cto, long int edges){   
    if(cfrom >= cto)
    {
        return;        
    }  
        
    if (*delta_q==NULL){
        node_delta_q *temp2 = (node_delta_q *)malloc(sizeof(node_delta_q));
        temp2->delta_q = val;
        temp2->from = cfrom;
        temp2->to = cto;        
        temp2->edges = edges; 
        temp2->left = temp2->right = NULL;        
        *delta_q = temp2;
        return;
    }
    
    if(val < (*delta_q)->delta_q){
        add_delta_q(&(*delta_q)->left, val, cfrom, cto, edges); 
    }  
    else{
        add_delta_q(&(*delta_q)->right, val, cfrom, cto, edges);
    }         
    return;    
}

void create_delta_q_row(node_delta_q **delta_q_row, node_edges_to *tree, long int from, const igraph_i_multilevel_community_list *communities){    
    if(tree==NULL)
    {   
        return;
    }    
    double delta_q = 2*(tree->edges - communities->item[from].weight_all*communities->item[tree->to].weight_all/(communities->weight_sum))/communities->weight_sum;
    
    // only create a node in the link list if delta_q>0
    if (delta_q>0) {
        add_delta_q(delta_q_row, delta_q, from, tree->to, tree->edges);        
    }    
    create_delta_q_row(delta_q_row, tree->left, from, communities);
    create_delta_q_row(delta_q_row, tree->right, from, communities);
}

// create a list of pairs of communities to be merged
void make_level_set(node_delta_q *tree, node_delta_q_list *level){    
    if(tree==NULL)
    {
        return;
    }   
    
    make_level_set(tree->right, level);    
    
    if ((tree->delta_q> 0) && (n_level<max_n_level)){
        level[n_level].delta_q = tree->delta_q;
        level[n_level].from = tree->from;
        level[n_level].to = tree->to;        
        level[n_level].edges = tree->edges;
        n_level++;
    }
    make_level_set(tree->left, level);
}

// convert a membership vector into communities
void membership_to_communities(const igraph_t *graph, const igraph_vector_t *weights, igraph_vector_t *membership, igraph_i_multilevel_community_list *communities){
  long int i;
  long int vcount = igraph_vcount(graph);
  long int ecount = igraph_ecount(graph);
  igraph_integer_t ffrom, fto;  
    
  communities->vertices_no = vcount;
  communities->communities_no = 0;
  communities->weights = weights;
  communities->weight_sum = 2 * igraph_vector_sum(weights);
  
  for (i=0; i < vcount; i++) {    
    communities->item[i].size = 0;
    communities->item[i].weight_inside = 0;
    communities->item[i].weight_all = 0;
  }
  
  for (i=0; i<vcount; i++){
      communities->item[(long int)VECTOR(*membership)[i]].size++;
      if (communities->item[(long int)VECTOR(*membership)[i]].size == 1) communities->communities_no++;
  }
  
  for (i = 0; i < ecount; i++) {
    igraph_real_t weight = 1;
    igraph_edge(graph, (igraph_integer_t) i, &ffrom, &fto);
    long int cfrom = (long int)VECTOR(*membership)[ffrom];
    long int cto = (long int)VECTOR(*membership)[fto];
        
    weight = VECTOR(*weights)[i]; 
    communities->item[cfrom].weight_all += weight;
    communities->item[cto].weight_all += weight;
    if (cfrom == cto)
      communities->item[cfrom].weight_inside += 2*weight;
  }

  communities->membership = membership;
  return;
}

void communities_backup(igraph_i_multilevel_community_list *communities_to, igraph_i_multilevel_community_list *communities_from, long int n){
    long int i;
    
    communities_to->vertices_no = communities_from->vertices_no;
    communities_to->communities_no = communities_from->communities_no;
    communities_to->weights = communities_from->weights;
    communities_to->weight_sum = communities_from->weight_sum;
    igraph_vector_update(communities_to->membership, communities_from->membership);
    for (i=0; i < n; i++) {      
        communities_to->item[i].size = communities_from->item[i].size;
        communities_to->item[i].weight_all = communities_from->item[i].weight_all;
        communities_to->item[i].weight_inside = communities_from->item[i].weight_inside;
    }  
}

int LPAm(const igraph_t *graph, igraph_vector_t *membership, igraph_i_multilevel_community_list *communities, igraph_vector_t *node_order){
  long int i, j;
  long int vcount = igraph_vcount(graph);
  igraph_real_t q, pass_q;
  int pass=0;
  igraph_bool_t changed = 0;
  
  igraph_vector_t links_community;
  igraph_vector_t links_weight;
  igraph_vector_t edges;
  
  igraph_vector_init(&links_community, 0);
  igraph_vector_init(&links_weight, 0);
  igraph_vector_init(&edges, 0);
  
  q = igraph_i_multilevel_community_modularity(communities);
  
  do { /* Pass begin */

      pass_q = q;
      changed = 0;
      long int index;
      for (index = 0; index < vcount; index++) {
          i = VECTOR(*node_order)[index];
                  
        /* Exclude vertex from its current community */
        igraph_real_t weight_all = 0;
        igraph_real_t weight_inside = 0;
        igraph_real_t weight_loop = 0;
        igraph_real_t max_q_gain = 0;
        igraph_real_t max_weight;
        long int old_id, new_id, n;

        igraph_i_multilevel_community_links(graph, communities, 
                                            (igraph_integer_t) i, &edges,
                                            &weight_all, &weight_inside, 
                                            &weight_loop, &links_community,
                                            &links_weight);

        old_id = (long int)VECTOR(*membership)[i];
        new_id = old_id;

        /* Update old community */
        igraph_vector_set(membership, i, -1);
        (*communities).item[old_id].size--;
        if (communities->item[old_id].size == 0) {communities->communities_no--;}
        communities->item[old_id].weight_all -= weight_all;
        communities->item[old_id].weight_inside -= 2*weight_inside + weight_loop;

         //debug("Remove %ld all: %lf Inside: %lf\n", i, -weight_all, -2*weight_inside + weight_loop); 

        /* Find new community to join with the best modification gain */
        max_q_gain = 0;
        max_weight = weight_inside;
        n = igraph_vector_size(&links_community);

        for (j = 0; j < n; j++) {
          long int c = (long int) VECTOR(links_community)[j];
          igraph_real_t w = VECTOR(links_weight)[j];

          igraph_real_t q_gain = 
            igraph_i_multilevel_community_modularity_gain(communities, 
                                                          (igraph_integer_t) c, 
                                                          (igraph_integer_t) i,
                                                          weight_all, w);
          
          /* debug("Link %ld -> %ld weight: %lf gain: %lf\n", i, c, (double) w, (double) q_gain); */
          if (q_gain > max_q_gain) {
            new_id = c;
            max_q_gain = q_gain;
            max_weight = w;
          }
        }

        /* debug("Added vertex %ld to community %ld (gain %lf).\n", i, new_id, (double) max_q_gain); */

        /* Add vertex to "new" community and update it */
        igraph_vector_set(membership, i, new_id);
        if (communities->item[new_id].size == 0) {communities->communities_no++;}
        communities->item[new_id].size++;
        communities->item[new_id].weight_all += weight_all;
        communities->item[new_id].weight_inside += 2*max_weight + weight_loop;

        if (new_id != old_id) {
          changed++;
        }
      }

      q = igraph_i_multilevel_community_modularity(communities);
      pass++;
      
    } while (changed && (q > pass_q)); /* Pass end */
  
  igraph_vector_destroy(&links_community);
  igraph_vector_destroy(&links_weight);
  igraph_vector_destroy(&edges);   
  
  IGRAPH_FINALLY_CLEAN(3);  
  
  return 0;
}

igraph_bool_t MetaLPAm_jump(const igraph_t *graph, igraph_i_multilevel_community_list *communities, igraph_real_t q_max, igraph_real_t delta, igraph_vector_t *node_order){
    long int i, j;
    long int vcount = igraph_vcount(graph);
    igraph_bool_t changed = 0;

    igraph_vector_t links_community;
    igraph_vector_t links_weight;
    igraph_vector_t edges;

    igraph_vector_t dominant_labels;
    igraph_vector_t dominant_weights;
    igraph_vector_t dominant_gains;
    
    igraph_vector_init(&links_community, 0);
    igraph_vector_init(&links_weight, 0);
    igraph_vector_init(&edges, 0);
    
    igraph_vector_init(&dominant_labels, 0);
    igraph_vector_init(&dominant_weights, 0);
    igraph_vector_init(&dominant_gains, 0);
    
    igraph_real_t q = igraph_i_multilevel_community_modularity(communities);    
    changed = 0;
    
    long int index;  
    for (index = 0; index < vcount; index++) {
        i = (long int)VECTOR(*node_order)[index];
      
        // get the list of neighbor communities of node i and weights inside its community
        igraph_real_t weight_all = 0;
        igraph_real_t weight_inside = 0;
        igraph_real_t weight_loop = 0;      

        igraph_i_multilevel_community_links(graph, communities, 
                                            (igraph_integer_t) i, &edges,
                                            &weight_all, &weight_inside, 
                                            &weight_loop, &links_community,
                                            &links_weight);
        
        long int old_id, new_id, n;         
        igraph_real_t q_gain=0;
        
        old_id = (long int)VECTOR(*(communities->membership))[i];
        new_id = old_id;         
       
        /* Find new community to join with following the Record-to-Record travel  */        
        igraph_vector_clear(&dominant_labels);
        igraph_vector_clear(&dominant_weights);
        igraph_vector_clear(&dominant_gains);
        
        n = igraph_vector_size(&links_community);        
        for (j = 0; j < n; j++) {
            long int c = (long int) VECTOR(links_community)[j];
            igraph_real_t w = VECTOR(links_weight)[j];            
            q_gain = 2*(w - weight_inside - weight_all/communities->weight_sum * (communities->item[c].weight_all - communities->item[old_id].weight_all + weight_all))/communities->weight_sum;
          
            if ((c!=old_id) && (( q + q_gain) > (q_max - delta))) {      
                igraph_vector_push_back(&dominant_labels, c);                
                igraph_vector_push_back(&dominant_weights, w); 
                igraph_vector_push_back(&dominant_gains, q_gain);   
            }
        }        
        
        if (igraph_vector_size(&dominant_labels) > 0) {
            
            /* Select randomly from the dominant labels */
            j = RNG_INTEGER(0, igraph_vector_size(&dominant_labels)-1);                       
            new_id = VECTOR(dominant_labels)[j];
            igraph_real_t weight = VECTOR(dominant_weights)[j];
            q_gain = VECTOR(dominant_gains)[j];
            
            /* Update community membership*/    
            igraph_vector_set(communities->membership, i, new_id);
            
            communities->item[old_id].size--;
            if (communities->item[old_id].size == 0) {communities->communities_no--;}
            communities->item[old_id].weight_all -= weight_all;
            communities->item[old_id].weight_inside -= 2*weight_inside + weight_loop;      
            
            if (communities->item[new_id].size == 0) {communities->communities_no++;}
            communities->item[new_id].size++;
            communities->item[new_id].weight_all += weight_all;
            communities->item[new_id].weight_inside += 2*weight + weight_loop;            
            
            q = q + q_gain;            
            changed++;           

            
        }     
        
      }
    
  igraph_vector_destroy(&links_community);
  igraph_vector_destroy(&links_weight);
  igraph_vector_destroy(&edges); 
  
  igraph_vector_destroy(&dominant_labels);
  igraph_vector_destroy(&dominant_weights);
  igraph_vector_destroy(&dominant_gains);
  
  IGRAPH_FINALLY_CLEAN(6);
  
  return changed;
}

int metaLPAm(const igraph_t *graph, igraph_i_multilevel_community_list *communities, igraph_real_t delta, int max_loop, igraph_vector_t *node_order){
  long int vcount = igraph_vcount(graph);
  int pass = 0;
  igraph_bool_t changed=0;
  
  igraph_vector_t membership_max;
  igraph_vector_init(&membership_max, vcount);      
  
  igraph_i_multilevel_community_list communities_max; // to return to the maximum found solution
  communities_max.item = igraph_Calloc(vcount, igraph_i_multilevel_community);
  communities_max.membership = &membership_max;
  communities_backup(&communities_max, communities, vcount);
  
  igraph_real_t q_max = igraph_i_multilevel_community_modularity(communities);
  int no_increase = 0;
  
  changed = MetaLPAm_jump(graph, communities, q_max, delta, node_order);
  pass++;
  
  while((changed) && (no_increase < max_loop)){
    pass = pass + LPAm(graph, communities->membership, communities, node_order);  
    no_increase++; 
    igraph_real_t q = igraph_i_multilevel_community_modularity(communities);
    if (q> q_max) {
        q_max = q;        
        communities_backup(&communities_max, communities, vcount);
        no_increase = 0;
    }         
    changed = MetaLPAm_jump(graph, communities, q_max, delta, node_order);
    pass++;    
  }   
  
  communities_backup(communities, &communities_max, vcount);
        
  igraph_vector_destroy(&membership_max);
  igraph_free(communities_max.item);
  
  return 0;
}

igraph_bool_t merging_communities(const igraph_t *graph, const igraph_vector_t *weights, igraph_i_multilevel_community_list *communities){    
    long int vcount = igraph_vcount(graph);
    long int ecount = igraph_ecount(graph);
    long int i, j;
    igraph_integer_t ffrom, fto;
            
    node_edges_to **edges_to;
    node_delta_q **delta_q_to;
    node_delta_q *delta_q_heap = NULL; 

    edges_to = (node_edges_to **) calloc (vcount, sizeof(node_edges_to*)); 
    delta_q_to = (node_delta_q **) calloc (vcount, sizeof(node_delta_q*)); 

    for (i=0; i< vcount; i++){
        edges_to[i] = NULL;
        delta_q_to[i] = NULL;
    };

    /* Edges between communities initialization :) */
    for (i = 0; i < ecount; i++) {
        igraph_real_t weight = 1;
        igraph_edge(graph, (igraph_integer_t) i, &ffrom, &fto);
        weight = VECTOR(*weights)[i]; 
        long int cfrom = (long int)VECTOR(*(communities->membership))[ffrom];
        long int cto = (long int)VECTOR(*(communities->membership))[fto];
        if (cfrom<cto) add_edges(&edges_to[cfrom], cto, weight);
        else if (cfrom>cto) add_edges(&edges_to[cto], cfrom, weight);
    }
  
    // create a binary tree of delta q for each row
    for (i=0; i< vcount; i++){
        create_delta_q_row(&delta_q_to[i], edges_to[i], i, communities);                
    }; 

    // creat a heap of maximal delta q
    for (i=0; i< vcount; i++){   
        node_delta_q *q = find_max(delta_q_to[i]);
        if (q!=NULL){
            add_delta_q(&delta_q_heap, q->delta_q, q->from, q->to, q->edges);
        }    
        free_edges_to(&edges_to[i]);
        free_delta_q(&delta_q_to[i]);        
    };

    // create a list of pairs of communities to merge
    node_delta_q_list *level;
    level = igraph_Calloc(max_n_level, node_delta_q_list);
    
    n_level = 0;  
    make_level_set(delta_q_heap, level);  
    free_delta_q(&delta_q_heap);
    igraph_bool_t touched[vcount]; 
    
    for(i=0; i<vcount; i++){
        touched[i] = 0;
    }
    
    for(i=0; i<n_level; i++){
        if ((touched[level[i].from] == 0) && (touched[level[i].to]==0) && (level[i].from!=level[i].to)){            
            for(j=0; j<vcount; j++){                      
                if ((long int)VECTOR(*(communities->membership))[j]== level[i].to){                        
                    VECTOR(*(communities->membership))[j] = level[i].from;                    
                }            
            }            
            touched[level[i].from]=1;
            touched[level[i].to]=1;
                    
            if (communities->item[level[i].from].size == 0) {communities->communities_no++;}
            communities->item[level[i].from].size += communities->item[level[i].to].size;
            communities->item[level[i].from].weight_all += communities->item[level[i].to].weight_all;
            communities->item[level[i].from].weight_inside += communities->item[level[i].to].weight_inside + 2*level[i].edges;            
            communities->item[level[i].to].size = 0;
            communities->communities_no--;
            communities->item[level[i].to].weight_all = 0;
            communities->item[level[i].to].weight_inside = 0;
        }
    }
    
    igraph_free(level);  
    IGRAPH_FINALLY_CLEAN(1);
    return n_level;    
}

int MetaLPAm_plus(const igraph_t *graph, const igraph_vector_t *weights, igraph_vector_t *membership, igraph_real_t *modularity) {    

    // data structure for the RRT inspired meta-heuristic
    const igraph_real_t delta = 0.01;
    const long int max_loop = 50;  

    long int vcount = igraph_vcount(graph);  
    igraph_bool_t changed=0; 
    igraph_vector_t w;

    if (weights) {
      IGRAPH_CHECK(igraph_vector_copy(&w, weights));   
      IGRAPH_FINALLY(igraph_vector_destroy, &w);  
    } else {
      IGRAPH_VECTOR_INIT_FINALLY(&w, igraph_ecount(graph));
      igraph_vector_fill(&w, 1);
    }

    igraph_vector_init_seq(membership, 0, vcount-1);  

    igraph_i_multilevel_community_list communities;
    communities.item = igraph_Calloc(vcount, igraph_i_multilevel_community);  
    membership_to_communities(graph, &w, membership, &communities); 

    // generate a random sequence for node order
    igraph_vector_t node_order;  
    igraph_vector_init_seq(&node_order, 0, vcount-1);
    igraph_vector_shuffle(&node_order);    

    LPAm(graph, membership, &communities, &node_order);    
    metaLPAm(graph, &communities, delta, max_loop, &node_order);  
    changed = merging_communities(graph, &w, &communities);   

    while(changed){          
        LPAm(graph, membership, &communities, &node_order);
        metaLPAm(graph, &communities, delta, max_loop, &node_order);  
        changed = merging_communities(graph, &w, &communities);  
    }     

    igraph_real_t q= igraph_i_multilevel_community_modularity(&communities);
  
    if (modularity) {    
      *modularity = q;
    }
  
  IGRAPH_CHECK(igraph_reindex_membership(membership, 0));
  
  igraph_free(communities.item); 
  igraph_vector_destroy(&node_order);
  
  IGRAPH_FINALLY_CLEAN(2);  
  return 0;
}

int main() {  
    int i;
    igraph_t g;    
    igraph_vector_t membership;
    igraph_real_t modularity;
    
    time_t t;    
    srand((unsigned) time(&t));
    
    igraph_vector_init(&membership, 0); 
    
    /* Zachary Karate club -- this is just a quick smoke test */
    igraph_small(&g, 0, IGRAPH_UNDIRECTED, 
                 0,  1,  0,  2,  0,  3,  0,  4,  0,  5,
                 0,  6,  0,  7,  0,  8,  0, 10,  0, 11,
                 0, 12,  0, 13,  0, 17,  0, 19,  0, 21,
                 0, 31,  1,  2,  1,  3,  1,  7,  1, 13,
                 1, 17,  1, 19,  1, 21,  1, 30,  2,  3,
                 2,  7,  2,  8,  2,  9,  2, 13,  2, 27,
                 2, 28,  2, 32,  3,  7,  3, 12,  3, 13,
                 4,  6,  4, 10,  5,  6,  5, 10,  5, 16,
                 6, 16,  8, 30,  8, 32,  8, 33,  9, 33,
                 13, 33, 14, 32, 14, 33, 15, 32, 15, 33,
                 18, 32, 18, 33, 19, 33, 20, 32, 20, 33,
                 22, 32, 22, 33, 23, 25, 23, 27, 23, 29,
                 23, 32, 23, 33, 24, 25, 24, 27, 24, 31,
                 25, 31, 26, 29, 26, 33, 27, 33, 28, 31,
                 28, 33, 29, 32, 29, 33, 30, 32, 30, 33,
                 31, 32, 31, 33, 32, 33,
                 -1);    
    
    MetaLPAm_plus(&g, 0, &membership, &modularity);
    
    printf("Modularity:  %.4f\n", modularity);
    printf("Membership: ");    
    for (i=0; i<igraph_vector_size(&membership); i++) {
         printf("%li ", (long int)VECTOR(membership)[i]);
    } 
            
    igraph_destroy(&g);    
    igraph_vector_destroy(&membership);  
    IGRAPH_FINALLY_CLEAN(2);  
    
    return 0;
}

