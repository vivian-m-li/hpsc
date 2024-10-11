To run:

- valgrind --tool=cachegrind --log-file="cachegrind.$1*$2.log" --cachegrind-out-file="cachegrind.$1*$2.out" ./main -n $1 -r $2

- cachegrind --tool=cachegrind --log-file="cachegrind.$1*$2.log" --cachegrind-out-file="cachegrind.$1*$2.out" ./main -n $1 -r $2

- kcachegrind cachegrind.$1\_$2.out

- All based on number of instructions that were performed, which is different than looking at the time
- Cache-misses can be involved so time is not necessarily correlated
- Look at Source File button in GUI --> this is how it's grouping the results (right next to the search query bar)
- You can expect the main function to be near 100%
- matVec was 68%
- sparse matrix format is much less full than the full matrix format -> can convert an array to sparse matrix format to
- improve performance
- arrayDouble was 31% (spent a lot of time allocating memory) IR refers to instruction reads so
- that's what we need to look at Relative button is depressed; shows absolute if you press it -> shows the actual count of the instructions
- Care about Data Read/Data Write Misses (to identify cache-misses), right of Shorten Templates
- Use LL (lower level) Data Read Miss
- The code is so bad the CPU can't even look ahead, one might be able to break the loop in standard matrix format into
- 10 different loops
- To run gprof, need to run the code with -p option for profiling -> output=gmon.out

Advice: take the code and make 1 change; keep the old capability next to the new capability for as long as possible with an if test that can pick either one, can gradually shrink different parts
