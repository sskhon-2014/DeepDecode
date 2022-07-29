![Alt text](logo.png?raw=true "Title")

## Description

This is a pure C++ instantiation of a decoder written for the primary purpose of decoding sequentially encoded molecular fluorescence data (eg. SeqFISH, SeqFISH+ data). This algorithm can, of course, be applied to diverse datasets such as astronomical observations. This is essentially a software that takes many observations/objects and clusters like objects. The primary purpose of this algorithm is to provide rapid, high-fidelity clustering of signal in 2D space ("decoding" of molecular idenities in SeqFISH data). The runtime of FastDecode in the analysis of observations in 2048x2048 space from sparse and high-density SeqFISH data reveals the runtime to be a few seconds to half-a-minute, rarely exceeding the 30 second mark.


## Test Run FastDecode

First, compile the C++ code:
```
g++ -o FastDecode *.cpp -std=c++11
```
Then, run FastDecode:
```
./FastDecode input.csv output.csv 45 81 4 2048 1
```
Result of clustering/decoding will appear as output.csv. See author for questions about customizing FastDecode to your desired application.


## Author

[@Harshaan Sekhon](https://www.linkedin.com/in/shaan-sekhon-1a217b154/)

## Version History

* 1.0, Initial Release

## License

Copyright (C) California Institute of Technology - All Rights Reserved
See LISCENCE.md for more details.
