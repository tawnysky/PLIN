Implementation of the paper "PLIN: A Persistent Learned Index for Non-Volatile Memory with High Read/Write Performance".

#### Dependence
We test our project on ubuntu LTS 2020. 
1. This project depends on libpmemobj from [PMDK](https://pmem.io/pmdk/libpmemobj/). Install it using command
    ```shell
    sudo apt-get install make
    sudo apt-get install libpmemobj-dev
    ``` 
2. Make sure you are avaliable with Optane memory (or you can use volatile memory to simulate pmem device on latest ubuntu version). The project assumes that your pmem device is mounted at address `/mnt/pmem` and you have write permission to it.

    (a). Avaiable to real Optane DC Memory, configure it in this [way](https://software.intel.com/content/www/us/en/develop/articles/qsg-part2-linux-provisioning-with-optane-pmem.html).
    
    (b). Not avaiable to Optane, try simluate it with DRAM following this [link](https://software.intel.com/content/www/us/en/develop/articles/how-to-emulate-persistent-memory-on-an-intel-architecture-server.html).


#### Usage
Follow these steps to play with PLIN.
1. Prepare the data file. Each line of the data file is a double type and contains at least 2e8 lines. We recommend using the normal_distribution from the STL library to generate the data. Configure the data file path in the benchmark function of test.cpp.
2. Configure the NVM pool file path in the benchmark function and the test function of test.cpp.
3. In the main function of test.cpp, choose to execute the test function or benchmark function.
4. Compile and execute. Use the command: "mkdir build; cd build; cmake .; make".

We presume that your machine is not capable of issuing clwb and clflushopt instrutions. If not, change the FLUSH_FLAG variable in makefile into -DCLWB and -DCLFLUSHOPT, respectively.
