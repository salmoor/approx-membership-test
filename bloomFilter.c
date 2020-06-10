/**
 * Author: Alemdar Salmoor
 * Date: December 15, 2019
 * Note: The hash functions are the implementations of the hash functions discussed in CLRS Book, Chapter 11.3 Hash Functions.
 * */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

struct Kmer
{

    //left ones correspond to 5' to 3'
    //right ones correspond to 3' to 5'
    unsigned long long left;
    unsigned long long right;
    unsigned long long h1left;
    unsigned long long h1right;
    unsigned long long h2left;
    unsigned long long h2right;
    unsigned long long h3left;
    unsigned long long h3right;
};




int NuctoDec(char c);
void readSequences(char * ref, char *** reference, int ** refSize, int * refCounter);
void createKmers( struct Kmer ** kmers, int * kmerCount, char ** reference, int * refSize, int numOfRefs, int kmerLen, long double maxValue, int slots, int prime, double aValue);
int smallerPrime(int x);
void setBit( int ** bloomFilter, unsigned long long k);
bool testBit( int * bloomFilter, unsigned long long k);
void populateBloom( int ** bloomFilter, struct Kmer * kmers, int kmerCount);
void testMembership(int * bloomFilter, struct Kmer * kmers, int kmerCount, int * matches);


int main(int argc, char** argv){

    char * ref = argv[2];
    char * queryF = argv[4];
    int kmerLen = atoi(argv[6]);
    int bloomsizeBytes = atoi(argv[8]);

    int bloomsizeBits = bloomsizeBytes * 8;

    char ** reference;
    int * refSize;
    int refCounter;

    //Read in references
    readSequences(ref, &reference, &refSize, &refCounter);


    //Read in Query
    char ** query;
    int * querySize;
    int queryCounter;

    readSequences(queryF, &query, &querySize, &queryCounter);

    double aValue = (sqrt(5) - 1)/2.0;
    // printf("A value: %lf\n", aValue);

    int smallerP = smallerPrime(bloomsizeBits);

    long double maxValue;
    maxValue = pow(4, kmerLen);

    
    struct Kmer * kmers;
    int kmerCount;
    createKmers(&kmers, &kmerCount, reference, refSize, refCounter, kmerLen, maxValue, bloomsizeBits, smallerP, aValue);

    struct Kmer * queryKmers;
    int queryKmerCount;

    createKmers(&queryKmers, &queryKmerCount, query, querySize, queryCounter, kmerLen, maxValue, bloomsizeBits, smallerP, aValue);

    
    double numberOfElems = (double) bloomsizeBits/32.0;
    int arrayLen = (int) ceil(numberOfElems);
    int * bloomFilter = malloc(arrayLen * sizeof(int));
    memset(bloomFilter, 0, arrayLen * sizeof(int));
    populateBloom(&bloomFilter, kmers, kmerCount);

    int matches = 0;

    testMembership(bloomFilter, queryKmers, queryKmerCount, &matches);

    printf("Number of k-mers indexed in reference: %d\n", kmerCount);
    printf("Number of k-mers scanned in query: %d\n", queryKmerCount);
    printf("Number of k-mers from query found in the reference: %d\n", matches);

    free(kmers);
    free(queryKmers);

    for (size_t i = 0; i < refCounter; i++)
    {
        free(reference[i]);
    }
    
    free(reference);
    free(refSize);

    for (size_t i = 0; i < queryCounter; i++)
    {
        free(query[i]);
    }

    free(query);
    free(querySize);

    free(bloomFilter);
    
    

    return 0;
}


void readSequences(char * ref, char *** reference, int ** refSize, int * refCounter){

    (* reference) = malloc(sizeof(char *));
    (* refSize) = malloc(sizeof(int));
    int * refCapacity = malloc(sizeof(int));

    int counterSize = 0;
    int counterCapacity = 1;

    int dsSize = 0;
    int dsCapacity = 1;

    //Read reference in
    (* refCounter) = 0;


    //input file for reference
    FILE * finputR = fopen(ref, "r");
    char letter;

    while (fscanf(finputR, "%c", &letter) > 0)
    {
        if (letter == '>')
        {
            name:
            while (fscanf(finputR, "%c", &letter) > 0)
            {
                if(letter == '\n'){
                    goto refs;
                }
            }
            
        }
        else
        {
            refs:
            if (dsSize >= dsCapacity)
            {
                (* reference) = realloc((* reference), 2 * dsCapacity * sizeof(char *));
                dsCapacity = dsCapacity * 2;
            }
        

            if(counterSize >= counterCapacity){
                (* refSize) = realloc((* refSize), 2 * counterCapacity * sizeof(int));
                refCapacity = realloc(refCapacity, 2 * counterCapacity * sizeof(int));
                counterCapacity = counterCapacity * 2;
            }

            (* reference)[(* refCounter)] = malloc(sizeof(char));
            (* refSize)[(* refCounter)] = 0;
            refCapacity[(* refCounter)] = 1;
            counterSize++;
            dsSize++;

            while (fscanf(finputR, "%c", &letter) > 0)
            {

                if (letter == '>'){
                    (* refCounter)++;
                    goto name;
                }

                if(letter != '\n'){
                    if ((* refSize)[(* refCounter)] >= refCapacity[(* refCounter)])
                    {
                        (* reference)[(* refCounter)] = realloc((* reference)[(* refCounter)], refCapacity[(* refCounter)] * 2 * sizeof(char));
                        refCapacity[(* refCounter)] = refCapacity[(* refCounter)] * 2;
                    }
                    
                    (* reference)[(* refCounter)][(* refSize)[(* refCounter)]] = letter;
                    (* refSize)[(* refCounter)]++;
                }

            }

        }
        
        
    }

    fclose(finputR);
    free(refCapacity);
    
}

void testMembership(int * bloomFilter, struct Kmer * kmers, int kmerCount, int * matches){

    int cnt = 0;

    bool response;

    while (cnt < kmerCount)
    {
        
        response = false;

        if (testBit( bloomFilter, kmers[cnt].h1left) && testBit( bloomFilter, kmers[cnt].h2left) && testBit( bloomFilter, kmers[cnt].h3left))
        {
            response = true;
            (* matches)++;
        }

        if (response == false)
        {
            if (testBit( bloomFilter, kmers[cnt].h1right) && testBit( bloomFilter, kmers[cnt].h2right) && testBit( bloomFilter, kmers[cnt].h3right))
            {
                (* matches)++;
            }
            
        }
        
        
        cnt++;
    }

}


void populateBloom( int ** bloomFilter, struct Kmer * kmers, int kmerCount){

    int cnt = 0;

    bool response;

    while (cnt < kmerCount)
    {
        
        response = false;

        if (testBit( (* bloomFilter), kmers[cnt].h1left) && testBit( (* bloomFilter), kmers[cnt].h2left) && testBit( (* bloomFilter), kmers[cnt].h3left))
        {
            response = true;
        }

        if (response == false)
        {
            if (testBit( (* bloomFilter), kmers[cnt].h1right) && testBit( (* bloomFilter), kmers[cnt].h2right) && testBit( (* bloomFilter), kmers[cnt].h3right))
            {
                response = true;
            }

            if (response == false)
            {
                setBit(bloomFilter, kmers[cnt].h1left);
                setBit(bloomFilter, kmers[cnt].h2left);
                setBit(bloomFilter, kmers[cnt].h3left);
                
            }
            
        }
        
        
        cnt++;
    }
    

}

void createKmers( struct Kmer ** kmers, int * kmerCount, char ** reference, int * refSize, int numOfRefs, int kmerLen, long double maxValue, int slots, int prime, double aValue){

    int kmerSize = 0;
    int kmerCapacity = 1;
    (* kmers) = malloc(sizeof(struct Kmer));

    int cnt = 0;
    int inner;
    int i;
    int j;
    unsigned long long sumLeft;
    unsigned long long sumRight;
    long double h1left;
    long double h1right;
    double diff;
    while (cnt <= numOfRefs)
    {
        if (refSize[cnt] >= kmerLen )
        {
            i = 0;
            while (i <= refSize[cnt] - kmerLen)
            {
                j = 0;
                sumLeft = 0;
                while (j < kmerLen)
                {

                    sumLeft = sumLeft + NuctoDec(reference[cnt][i + j]) * (unsigned long) pow(4, j);
                    j++;
                }
                j--;

                inner = 0;
                sumRight = 0;
                while (j >= 0)
                {
                    sumRight = sumRight + NuctoDec(reference[cnt][i + j]) * (unsigned long) pow(4, inner);
                    inner++;
                    j--;
                }

                if(kmerSize >= kmerCapacity){
                    (* kmers) = realloc((* kmers), 2 * kmerCapacity * sizeof(struct Kmer));
                    kmerCapacity = kmerCapacity * 2;
                }
                (* kmers)[kmerSize].left = sumLeft;
                (* kmers)[kmerSize].right = sumRight;
                h1left = (double) (sumLeft * 1.0)/maxValue;
                h1right = (double) (sumRight * 1.0)/maxValue;

                (* kmers)[kmerSize].h1left = (unsigned long long) floor(h1left * slots);
                (* kmers)[kmerSize].h1right = (unsigned long long) floor(h1right * slots);

                (* kmers)[kmerSize].h2left = sumLeft % prime;
                (* kmers)[kmerSize].h2right = sumRight % prime;

                diff = sumLeft * aValue - floor(sumLeft * aValue);

                //printf("sumLeft * aValue: %lf, diff: %lf\n", sumLeft * aValue, diff);
                (* kmers)[kmerSize].h3left = (unsigned long long) floor(slots * diff);

                diff = sumRight * aValue - floor(sumRight * aValue);

                //printf("sumRight * aValue: %lf, diff: %lf\n", sumRight * aValue, diff);
                (* kmers)[kmerSize].h3right = (unsigned long long) floor(slots * diff);

                
                kmerSize++;


                i++;
                
            }
            
        }

        cnt++;
        
    }

    *kmerCount = kmerSize;
    

}

int smallerPrime(int x){

    int i, n;
    if( x % 2  == 0){ n = x - 1; }
    else{ n = x - 2; }

    bool prime = true;

    while (true)
    {
        i = 2;
        while (i < n/2)
        {
            if (n%i == 0) { prime = false; break; }

            i++;
        }
        if (prime) { break; }
        n = n - 2;
        prime = true;
        
    }

    return n;
    
}

void setBit( int ** bloomFilter, unsigned long long k){

    unsigned long long arrIndex = k/32;

    unsigned long long index = k%32;

    unsigned int flag = 1;

    flag = flag << index;

    (*bloomFilter)[arrIndex] = (*bloomFilter)[arrIndex] | flag;

}

bool testBit( int * bloomFilter, unsigned long long k){

    unsigned long long arrIndex = k/32;

    unsigned long long index = k%32;

    unsigned int flag = 1;

    flag = flag << index;

    if ( bloomFilter[arrIndex] & flag ){
        return true;
    }else{
        return false;
    }

}

int NuctoDec(char c){
    if (c == 'A'){ return 0; }
    else if(c == 'C'){ return 1; }
    else if(c == 'G'){ return 2; }
    else { return 3; }
}