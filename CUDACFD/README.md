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
