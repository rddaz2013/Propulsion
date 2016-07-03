## CUDA explained

cf. [CUDA: cudaEvent_t and cudaThreadSynchronize usage](http://stackoverflow.com/questions/5801717/cuda-cudaevent-t-and-cudathreadsynchronize-usage)

```
cudaEvent_t start, stop;

cudaEventRecord( start, 0 ); // where 0 is the default stream

convect<<<grid,block>>>(...) // also using the default stream

cudaEventRecord( stop,  0 ); // also using the default stream

// synchronize
cudaEventSynchronize(start); // optional
cudaEventSynchronize(stop) ; // wait for the event to be executed completely!

// calculate elapsed time
float dt_ms;
cudaEventElapsedTime( &dt_ms, start, stop);

```

##### 20160703

[GCC compile error with >2 GB of code](http://stackoverflow.com/questions/6296837/gcc-compile-error-with-2-gb-of-code)

This question applied to my case because I was obtaining errors in compiling such as the following:

```  
/tmp/ccbxDNJI.o: In function `__static_initialization_and_destruction_0(int, int)':
convect3dupwind_shared.cpp:(.text+0x52): relocation truncated to fit: R_X86_64_PC32 against `.bss'
collect2: error: ld returned 1 exit status  
```  

In particular, **relocation truncated to fit**.  

[BЈовић](http://stackoverflow.com/users/476681/b%d0%88%d0%be%d0%b2%d0%b8%d1%9b) had an answer that seemed to apply:

> The [x86-64 ABI](http://www.x86-64.org/documentation/abi.pdf) used by Linux defines a "large model" specifically to avoid such size limitations, which includes 64-bit relocation types for the GOT and PLT. (See the table in section 4.4.2, and the instruction sequences in 3.5.5 which show how they are used.)

> Since your functions are occupying 2.8 GB, you are out of luck, because gcc doesn't support large models. What you can do, is to reorganize your code in such a way that would allow you to split it into shared libraries which you would dynamically link.

> If that is not possible, as someone suggested, instead of putting your data into code (compiling and linking it), since it is huge, you can load it at run time (either as a normal file, or you can mmap it).

> EDIT

> Seems like the large model is supported by gcc 4.6 (see this page). You can try that, but the above still applies about reorganizing your code.

This is because when I "dialed down" the "resolution" or grid size from $N_x \times N_y \times N_z = 1920 \times 1920 \times 32 = 117964800$  I dialed it down to $N_x \times N_y \times N_z = 1080 \times 1080 \times 64 = 74649600$

This note is for file, pertains to, `convect3dupwind_shared.cpp`



