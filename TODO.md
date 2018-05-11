### TODO list 

- Add output file reporting the options that were used to run the code, e.g. _number of significant clusters_, etc.

- Replace warning functions from class `SPEDisplay` with a unique
  function `msg = f(mainfilename)`.

- The `bang` character (`!`) is a shell-escape only valid for UNIX-like system. While trying to call the external executable `graphMetrics` in mainDRTGraphMetrics in Windows, an error occurs. Not sure if due to POSIX or compilation. 

- Code requires an entire review to improve performance, remove cluttering and bad programming. Taking too long to finish for SPE graphMetrics. 

- Make an output report file to log the execution parameters. 