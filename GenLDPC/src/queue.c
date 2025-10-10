#include "queue.h"
#include <stdlib.h>
#include <stdio.h>

#define MinQueueSize (5)

struct QueueRecord
{
    int Capacity;
    int Front;
    int Rear;
    int Size;
    ElementType *Array;
};

/* Check if the queue is empty */
int IsEmpty(Queue Q)
{
    return Q->Size == 0;
}

/* Check if the queue is full */
int IsFull(Queue Q)
{
    return Q->Size == Q->Capacity;
}

/* Create a queue with a maximum size */
Queue CreateQueue(int MaxElements)
{
    Queue Q;

    if (MaxElements < MinQueueSize)
    {
        fprintf(stderr, "Queue size is too small\n");
        return NULL;
    }

    Q = (Queue)malloc(sizeof(struct QueueRecord));
    if (Q == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        return NULL;
    }

    Q->Array = (ElementType *)malloc(sizeof(ElementType) * MaxElements);
    if (Q->Array == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        free(Q);
        return NULL;
    }

    Q->Capacity = MaxElements;
    MakeEmpty(Q);

    return Q;
}

/* Make the queue empty */
void MakeEmpty(Queue Q)
{
    Q->Size = 0;
    Q->Front = 1;
    Q->Rear = 0;
}

/* Dispose the queue */
void DisposeQueue(Queue Q)
{
    if (Q != NULL)
    {
        free(Q->Array);
        free(Q);
    }
}

static int Succ(int Value, Queue Q)
{
    if (++Value == Q->Capacity)
        Value = 0;
    return Value;
}

/* Add an element to the queue */
void Enqueue(ElementType X, Queue Q)
{
    if (IsFull(Q))
    {
        fprintf(stderr, "Full queue\n");
        return;
    }
    else
    {
        Q->Size++;
        Q->Rear = Succ(Q->Rear, Q);
        Q->Array[Q->Rear] = X;
    }
}

/* Get the front element of the queue */
ElementType Front(Queue Q)
{
    if (IsEmpty(Q))
    {
        fprintf(stderr, "Empty queue\n");
        return 0;
    }
    return Q->Array[Q->Front];
}

/* Remove the front element from the queue */
void Dequeue(Queue Q)
{
    if (IsEmpty(Q))
    {
        fprintf(stderr, "Empty queue\n");
        return;
    }
    else {
        Q->Size--;
        Q->Front = Succ(Q->Front, Q);
    }
}

/* Get and remove the front element from the queue */
ElementType FrontAndDequeue(Queue Q)
{
    ElementType X = 0;

    if (IsEmpty(Q))
    {
        fprintf(stderr, "Empty queue\n");
        return 0;
    }

    Q->Size--;
    X = Q->Array[Q->Front];
    Q->Front = Succ(Q->Front, Q);

    return X;
}
