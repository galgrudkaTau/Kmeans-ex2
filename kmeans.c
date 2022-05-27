#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

typedef struct list ELEMENT;
typedef ELEMENT *LINK;
typedef double *datapoint; 
/* LINK is a pointer to ELEMENT */
/* datapoint is a pointer to a double array */

static void anErrorHasOccurred();
static void invalidInput();
static void restartClusters(LINK *clusters, int k, int conti);
static void delete_list(LINK head);
static void kMeans(int k, int size, int d,double epsilon, int max_iter, double **cents, LINK *clusters, double **matrix);
static void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d);
static int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d, double epsilon);
static double calculateNorma(double *old, double * new, int d);
static double calculateDistance(double *datapoint, double *centroid, int d);
static double** kMeansMain(int size, int k, int d ,int max_iter,double epsilon ,PyObject *cents,PyObject *datapoints);
static PyObject* fit(PyObject *self,PyObject *args);

struct list { 
    int datapoint; 
    struct list *next;   
};

static PyMethodDef KmeansCAPIMethods[]={
    {"fit",
     (PyCFunction)fit,
     METH_VARARGS,
     PyDoc_STR("")},
     {NULL, NULL,0,NULL}
    };

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp", 
        NULL, 
        -1, 
        KmeansCAPIMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void){
    PyObject *m;
    m=PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
    }

static PyObject* fit(PyObject *self,PyObject *args){
    printf("entered fit\n");
    int max_iter, k, size, d;
    double epsilon;
    PyObject* initCentArray;
    PyObject* inputMatrix;
    if(!PyArg_ParseTuple(args, "OOiiiid",&initCentArray, &inputMatrix, &size, &k, &d ,&max_iter, &epsilon)){
        return NULL;
    }
     
    return Py_BuildValue("O",kMeansMain(size, k, d ,max_iter,epsilon,initCentArray,inputMatrix));
}

static void anErrorHasOccurred(){
    printf("An Error Has Occurred\n");
    exit(1);
}

void invalidInput(){
    printf("Invalid Input!\n");
    exit(1);
}


static void restartClusters(LINK *clusters, int k, int conti) {
    int i;
    LINK current;
    if (conti>0){ /* should do restart and new intilazation */
        for (i=0; i<k; i++) {
        current = clusters[i];
        delete_list(current);
        clusters[i] = (ELEMENT*)malloc(sizeof(ELEMENT));
        if (clusters[i] == NULL){
            anErrorHasOccurred();
        }
        clusters[i]->datapoint = -1;
        clusters[i]->next = NULL;
        }
    }
    else{
       for (i=0; i<k; i++) {
        current = clusters[i];
        delete_list(current);
        } 
    } 
}

static void delete_list(LINK head) {
    if (head != NULL){
        delete_list(head->next);
        free(head);
    }
}

static void kMeans(int k, int size, int d, double epsilon, int max_iter, double **cents, LINK *clusters, double **matrix){
    int iter = 0;
    int continue_condition = 1;
    while (iter<max_iter && continue_condition){
       assignToCluster(cents, matrix, clusters, k, size, d);
       continue_condition = updateCentroids(cents, clusters, matrix, k, d, epsilon);
       iter++;
    }
}

static void assignToCluster(double **cents,double** datapoints,LINK *clusters, int k, int size, int d) {
    double distance;
    double tempDist = 0;
    int clusterIndex, i, j;
    LINK current, new;
 
    for(i=0; i<size; i++) {
        clusterIndex = -1;
        distance = DBL_MAX;
        for (j=0; j<k; j++){
            tempDist = calculateDistance(datapoints[i], cents[j], d);
            if (tempDist < distance) {
                distance = tempDist;
                clusterIndex = j;
            }
        }
        if (clusters[clusterIndex]->datapoint == -1) {
            clusters[clusterIndex]->datapoint = i;
        }else {
            current = clusters[clusterIndex];
            new = (ELEMENT*)malloc( sizeof(ELEMENT));
            if (new == NULL){
                anErrorHasOccurred();
            }
            new -> datapoint = i;
            new -> next = current;
            clusters[clusterIndex] = new;
        }
    }
}

static int updateCentroids(double **cents, LINK *clusters, double** inputMatrix ,int k, int d, double epsilon) {
    int i, j, m, s, difference = 0, sizeOfCluster;
    LINK current = NULL;
    double *sum = (double *)calloc(d, sizeof(double));
    if (sum == NULL){
        anErrorHasOccurred();
    }
    for (i=0; i<k; i++) {
        current = clusters[i];
        sizeOfCluster = 0;
        while (current != NULL && current->datapoint != -1){
            for  (j=0; j<d; j++) {
                sum[j] += inputMatrix[current->datapoint][j];
            }
            current = current->next;
            sizeOfCluster++;
        }
        if (sizeOfCluster == 0){
            anErrorHasOccurred();
        }
        for (m=0; m<d; m++){
            sum[m] /= sizeOfCluster;
        }
        difference += (calculateNorma(cents[i], sum, d) > epsilon)?1:0;
        for (s=0; s<d; s++) {
            cents[i][s] = sum[s];
            sum[s] = 0;
        }
    }
    restartClusters(clusters, k, difference);
    free(sum);
    return difference > 0;
}

static double calculateNorma(double *old, double *new, int d){
    double sum = 0;
    int i;
    for (i=0; i<d; i++){
        sum += pow(old[i]-new[i],2);
    }
    return sqrt(sum);
}

static double calculateDistance(double *datapoint, double *centroid, int d){
    double distance = 0;
    int j;
    for (j = 0; j<d; j++){
        distance+= pow(datapoint[j] - centroid[j],2);
    }
    return distance;
}

static double **objectToMatrix(PyObject* obj, int amount, int length){
    int i,j;
    double *vector = NULL;
    double **matrix = NULL;    PyObject *item = NULL;
    printf("items in array =%d\n", amount);
    printf("length of item =%d\n", length);
    vector = (double *)calloc(amount*length, sizeof(double));
    if (vector == NULL){
        printf("error in vector allocation\n");
        anErrorHasOccurred();
    }
    matrix = (double **)calloc(length, sizeof(double *)); 
    if (matrix == NULL){
        printf("error in vector allocation\n");
        anErrorHasOccurred();
    }
    for (i=0; i<amount; i++) {
        matrix[i] = vector + i*length;
    }
    matrix[0][0]= 1.0;
    printf("about to enter for loop\n");
    for (i=0 ; i<amount ; i++){
        item = PyList_GetItem(obj,i);
        for (j=0; j<length; j++){

            matrix[i][j]=PyFloat_AsDouble(PyList_GetItem(item,j));
        }
    }
    printf("reached end of Object to Matrix method\n");
    return matrix;
}

static double **kMeansMain(int size, int k, int d, int max_iter, double epsilon ,PyObject* cents,PyObject* datapoints){
    int i, j;
    double ** centroids, **dataMatrix;
    LINK *clusters;
    /*convert PyObject to list,..*/
    printf("entered C main\n");
    centroids = objectToMatrix(cents, k, d);
    dataMatrix = objectToMatrix(datapoints, size, d);
    /*d = sizeof(*centroids)[0]/ sizeof(double);
    size = sizeof(*dataMatrix)/sizeof(*centroids)[0];
    k = sizeof(*centroids)/sizeof(*centroids)[0];*/
    clusters = (LINK *)calloc(k, sizeof(LINK));
    if(clusters == NULL){
        anErrorHasOccurred();
    }
    printf("allocated clusters\n");
    restartClusters(clusters, k, 1); /*this array holds K datapoints*/
    printf("initialized clusters\n");
    kMeans(k, size, d, epsilon, max_iter, centroids, clusters, dataMatrix);
    printf("ran kMeans\n");
    free(dataMatrix[0]); 
    free(dataMatrix);
    free(clusters); 
    
    for (i = 0; i<k; i++){
        for (j=0; j<d-1; j++){
            printf("%.4f,",centroids[i][j]);
        }
        printf("%.4f",centroids[i][j]);
        printf("\n");
    }
    printf("finished code");
    return centroids;
}
