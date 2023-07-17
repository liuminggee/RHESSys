#include <stdlib.h>
#include "rhessys.h"
#include <stdio.h>
#include <float.h>
#include <math.h>
/*------------------------------
sort heights in descending order
-------------------------------*/  
int key_compare( void * e1,  void *e2 )
{
	/*------------------------------------------------------*/
	/*	Local Function Definition. 							*/
	/*------------------------------------------------------*/
	double v1,v2;	
	/*------------------------------------------------------*/
	/*	Local Variable Definition. 							*/
	/*------------------------------------------------------*/
		v1 = ((struct layer_object *)e1)->height;
		v2 = ((struct layer_object *)e2)->height;

	return (v1>v2) ? -1 : (v1<v2) ? 1: 0;
}/*end key_compare.c*/
//05202022LML compare double
//int close_enough(double a, double b)
//{
//    if (fabs(a - b) <= DBL_EPSILON * fmax(fabs(a), fabs(b))) return 1;
//    return 0;
//}

//For calculating moving average
void initialize_FIFO_Queue(FIFO_Queue* queue, int maxisize)
{
  queue->front = 0;
  queue->rear = 0;
  queue->size = 0;
  queue->maxsize = maxisize;
}
//______________________________________________________________________________
int isEmpty_FIFO_Queue(FIFO_Queue* queue)
{
  return queue->front == 0;
}
//______________________________________________________________________________
double dequeue_FIFO_Queue(FIFO_Queue* queue)
{
  if (isEmpty_FIFO_Queue(queue)) {
    printf("Queue is empty. Cannot dequeue element.\n");
    return -1;
  }

  double value = queue->front->data;
  FIFO_Node* temp = queue->front;
  queue->front = queue->front->next;
  queue->size--;
  free(temp);
  if (queue->front == 0) {
    queue->rear = 0; // Queue is now empty
    queue->size = 0;
  }
  return value;
}
//______________________________________________________________________________
void enqueue_FIFO_Queue(FIFO_Queue* queue, double value)
{
  FIFO_Node* newNode = (FIFO_Node*) malloc(sizeof(FIFO_Node));
  newNode->data = value;
  newNode->next = 0;
  if (isEmpty_FIFO_Queue(queue)) {
    queue->front = newNode;
    queue->rear = newNode;
  } else {
    queue->rear->next = newNode;
    queue->rear = newNode;
  }
  queue->size++;
  if (queue->size > queue->maxsize)
      dequeue_FIFO_Queue(queue);
}
//______________________________________________________________________________
int peek_FIFO_Queue(FIFO_Queue* queue)
{
  if (isEmpty_FIFO_Queue(queue)) {
    printf("Queue is empty. Cannot peek element.\n");
    return -9999;
  }
  return queue->front->data;
}
//______________________________________________________________________________
void printQueue_FIFO_Queue(FIFO_Queue* queue)
{
  if (isEmpty_FIFO_Queue(queue)) {
    printf("Queue is empty.\n");
    return;
  }
  printf("Queue: ");
  FIFO_Node* current = queue->front;
  while (current != 0) {
    printf("%f ", current->data);
    current = current->next;
  }
  printf("\n");
}
//______________________________________________________________________________
double avgvalue_FIFO_Queue(FIFO_Queue* queue)
{
  if (isEmpty_FIFO_Queue(queue)) {
    //printf("Queue is empty.\n");
    return 0;
  }
  FIFO_Node* current = queue->front;
  double avg = 0;
  while (current != 0) {
    avg += current->data;
    current = current->next;
  }
  return avg/(double)queue->size;
}
