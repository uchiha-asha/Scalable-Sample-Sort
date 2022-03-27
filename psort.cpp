#include "psort.h"
#include <omp.h>
#include <iostream>
#include <string.h>
#include <chrono>

// scratch/pbs/pbs.2628575.pbshpc.x8z

const uint32_t bits = 8, base = (1<<bits)-1;

void InsertionSort(uint32_t *data, uint32_t n)
{
    for (uint32_t i=1; i<n; i++) {
        int j = i-1;
        uint32_t key=data[i];
        while (j >= 0 && data[j] > key) {
            data[j+1] = data[j];
            j--;
        }
        data[j+1] = key;
    }
}

void CountSort(uint32_t *data, uint32_t n, uint32_t exp)
{
    uint32_t *output = new uint32_t[n];
    uint32_t *count = new uint32_t[(1<<bits)];

    for (uint32_t i=0; i<=base; i++)
        count[i] = 0;

    for (uint32_t i=0; i<n; i++) {
        count[(data[i] >> exp) & base]++;
    }

    for (uint32_t i=1; i<=base; i++) {
        count[i] += count[i-1];
    }

    for (int64_t i=n-1; i>=0; i--) {
        // std::cout << count[(data[i] >> exp) & base] << " " << ((data[i] >> exp) & base) << std::endl;
        output[count[(data[i] >> exp) & base] - 1] = data[i];
        count[(data[i] >> exp) & base]--;
    }
    // std::cout << "here" << std::endl;
    for (uint32_t i=0; i<n; i++) {
        data[i] = output[i];
    }

    delete[] output;
    delete[] count;
}

void SequentialSort(uint32_t *data, uint32_t n)
{
    // std::sort(data, data+n);
    // std::cout << "Sequential Sorting: " << n << std::endl;
    if (n <= 16) {
        InsertionSort(data, n);
        return;
    }

    for (uint32_t i=0; i<32; i+=bits) {
        CountSort(data, n, i);
    }
}

void getSplitters(uint32_t *data, uint32_t *splitters, uint32_t n, int p)
{
    uint32_t *pseudo_splitters = new uint32_t[p*p];
    for (int i=0; i<p; i++) {
        for (int j=0; j<p; j++) {
            pseudo_splitters[i*p+j] = data[i*(n/p) + j];
        }
    }
    SequentialSort(pseudo_splitters, p*p);
    splitters[0] = 0, splitters[p] = UINT32_MAX;
    for (int i=1; i<p; i++) {
        splitters[i] = pseudo_splitters[i*p];
    }
    /*for (int i=1; i<=p; i++) std::cout << splitters[i] << " ";
    std::cout << std::endl;*/
    delete[] pseudo_splitters;
}

void create_bst(uint32_t *bst, uint32_t *splitters, int index, int l, int r)
{
    if (l > r) 
        return;
    uint32_t mid = (l + r)/2;
    // std::cout << mid << " " << l << " " << r << " " << index << " " << splitters[mid] << std::endl;
    bst[index] = splitters[mid];
    create_bst(bst, splitters, 2*index, l, mid-1);
    create_bst(bst, splitters, 2*index+1, mid+1, r);
}

inline uint8_t get_bucket(uint32_t &val, uint32_t *splitters, int &p)
{
    uint8_t bucket;
    // bucket = std::upper_bound(splitters, splitters+p+1, val) - splitters;
    for (bucket=1; bucket<=p; bucket++) {
        if (val > splitters[bucket-1] && val <= splitters[bucket]) {
            return bucket;
        }
    }
    return bucket - (bucket==p+1);
}

inline uint8_t get_bucket_optimized(uint32_t &val, uint32_t *bst, int &sz, int &bit_sz)
{
    int j=1;
    #pragma GCC unroll 6
    for (int i=0; i<bit_sz; i++) {
        // std::cout << "here" << val << " " << j << " " << i << " " << bit_sz << " " << sz << std::endl;
        j = 2*j + (val > bst[j]);
    }
    // std::cout << "here" << j-sz+1 << std::endl;
    return j-sz+1;
}

void assign_to_bucket(uint32_t *data, uint32_t **bucket_n_thread, uint32_t *splitters, 
                        uint8_t *assigned_buckets, uint32_t &n, int &p,
                        uint32_t *bst, int &sz, int &bit_sz, int tid)
{
    if (p <= 8) {
        for (uint64_t i=0; i<n; i++) {
            uint8_t bucket = get_bucket(data[i], splitters, p);
            assigned_buckets[i] = bucket;
            // #pragma omp atomic
            bucket_n_thread[tid][bucket]++;
        }
    }
    else {
        for (uint64_t i=0; i<n; i++) {
            uint8_t bucket = get_bucket_optimized(data[i], bst, sz, bit_sz);
            assigned_buckets[i] = bucket;
            // #pragma omp atomic
            bucket_n_thread[tid][bucket]++;
        }
    }
    
}
// int l=0;
void ParallelSort(uint32_t *data, uint32_t n, int p)
{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.
    // if (l++ > 5) return;
    // std::cout << "Parallel Sorting: " << n << " " << p << std::endl;
    if ((uint64_t)p*p >= (uint64_t)2*n || n < (1<<6)) {
        SequentialSort(data, n);
        return;
    }
    
    uint32_t *splitters = new uint32_t[2*p+1];
    uint32_t *data_in_buckets = new uint32_t[n];
    uint8_t *assigned_buckets = new uint8_t[n];
    uint32_t *bucket_n = new uint32_t[p+1];
    uint32_t *bst = new uint32_t[2*p+1];
    uint32_t *global_bucket_start = new uint32_t[p+1];

    int n_threads = omp_get_num_threads();

    uint32_t **bucket_n_thread = new uint32_t*[n_threads];
    uint32_t **bucket_start = new uint32_t*[n_threads];
    splitters[0] = 0;
    for (int i=0; i<=p; i++) {
        bucket_n[i] = 0;
        splitters[p+i] = UINT32_MAX;
    }

    for (int i=0; i<n_threads; i++) {
        bucket_n_thread[i] = new uint32_t[p+1];
        bucket_start[i] = new uint32_t[p+1];
        for (int j=0; j<=p; j++) {
            bucket_n_thread[i][j] = 0;
        }
    }

    int sz = 1, bit_sz = 0;
    while (sz < p) sz = (sz << 1), bit_sz++;

    getSplitters(data, splitters, n, p);

    create_bst(bst, splitters+1, 1, 1, sz-1);

    for (int i=0; i<n_threads; i++) {
        #pragma omp task
        {
            uint32_t bucket_size = n/n_threads, start = bucket_size*i;
            if (i==n_threads-1) {
                bucket_size = n - bucket_size*i;
            }

            assign_to_bucket(data+start, bucket_n_thread, splitters, assigned_buckets+start, bucket_size, p, bst, sz, bit_sz, i);
            uint32_t *added_to_buckets = new uint32_t[p+1];
            added_to_buckets[1] = 0, bucket_start[i][1] = start;
            for (int j=2; j<=p; j++) {
                bucket_start[i][j] = bucket_start[i][j-1] + bucket_n_thread[i][j-1];
                added_to_buckets[j] = 0;
            }
            for (uint64_t j=0; j<bucket_size; j++) {
                data_in_buckets[bucket_start[i][assigned_buckets[start+j]] + (added_to_buckets[assigned_buckets[start+j]]++)] = data[start+j];
            }
        }
    }

    #pragma omp taskwait
    
    for (int i=0; i<n_threads; i++) {
        for (int j=0; j<=p; j++) {
            bucket_n[j] += bucket_n_thread[i][j];
        }
    }
    global_bucket_start[1] = 0;
    for (int i=2; i<=p; i++) 
        global_bucket_start[i] = global_bucket_start[i-1] + bucket_n[i-1];

    // for (int i=1; i<=p; i++) std::cout << global_bucket_start[i] << " " << bucket_n[i] << std::endl;
    for (int i=1; i<=p; i++) {
        #pragma omp task
        {
            uint32_t curt = 0, count = 0;
            for (uint64_t j=0; j<bucket_n[i]; j++) {
                while (curt < (uint32_t)n_threads && count==bucket_n_thread[curt][i]) curt++, count=0;
                if (curt == (uint32_t)n_threads) break;
                data[global_bucket_start[i]+j] = data_in_buckets[bucket_start[curt][i]+(count++)];
            }
            if ((uint64_t)bucket_n[i] < (uint64_t)2*n/p)
                SequentialSort(data+global_bucket_start[i], bucket_n[i]);
            else 
                ParallelSort(data+global_bucket_start[i], bucket_n[i], p);
        }
    }

    #pragma omp taskwait

    delete[] splitters;
    delete[] data_in_buckets;
    delete[] assigned_buckets;
    /*std::cout << "ha" << std::endl;
    for (int i=1; i<=p; i++) {
        #pragma omp task
        {
            if ((uint64_t)bucket_n[i] < (uint64_t)2*n/p)
                SequentialSort(data+global_bucket_start[i], bucket_n[i]);
            else 
                ParallelSort(data+global_bucket_start[i], bucket_n[i], p);
        }
    }

    #pragma omp taskwait*/

    delete[] global_bucket_start;
    delete[] bucket_start;
    delete[] bucket_n;
}
