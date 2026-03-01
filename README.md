# Lower-Bound Distance Queries under Partial Information (VLDB 2026)
Code repository for metric space lower bound query paper submitter in VLDB 2026

Please make sure to have boost c++ library and openmp support librares installed in the system. 

Please compile the code using the following command: 

```
g++ -std=c++23 -O3 -DNDEBUG -march=native -pthread *.cpp -o edge_landmarks_sampling
```

For execution of a simple auto created uniform instance graph please run the following command. 

```
./edge_landmarks_sampling
```

This will execute the set of experiments assuming the default settings provided in the code. 
