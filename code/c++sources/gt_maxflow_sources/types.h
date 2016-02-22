#ifndef hi_pr_types_h
#define hi_pr_types_h

typedef  /* arc */
   struct arcSt
{
    float           resCap;          /* residual capacity */
    struct nodeSt   *head;           /* arc head */
    struct arcSt    *rev;            /* reverse arc */
}
    arc;

typedef  /* node */
   struct nodeSt
{
    arc *first;            /* first outgoing arc */
    arc *current;          /* current outgoing arc */
    double excess;         /* excess at the node */
    long d;                /* distance label */
    struct nodeSt *bNext;  /* next node in bucket */
    struct nodeSt *bPrev;  /* previous node in bucket */
} node;


typedef /* bucket */
   struct bucketSt
{
  node *firstActive;      /* first node with positive excess */
  node *firstInactive;    /* first node with zero excess */
} bucket;

#endif
