![Alt text](logo.png?raw=true "Title")

## Description

DeepDecode is a deep-learning algorithm written for the primary purpose of decoding sequentially encoded molecular fluorescence data (eg. SeqFISH, SeqFISH+ data). Given several 2D images that represent dispirate "barcoding rounds," this algorithm can essentially detect signal and determine parity all in one step. This takes much of the guesswork out of signal detection and decoding of SeqFISH datasets.

## Neural Network Architecture
![Alt text](model.png)

## DeepDecoded Dataset Result - Mouse Testes
### High-Resolution Tissue Segmentation
This is a completely decoded SeqFISH dataset - the mouse testes - generated using DeepDecode.

![Alt text](mouse-testes-segmentation.png)

### Clustering and Gene Enrichment Analysis
The following are the results of UMAP generation and clustering applied to the aforementioned dataset, which produced the image you see above:

![Alt text](mouse-testes-parameters.png)

### Credit: Experiment performed by the Cai Lab at Caltech.

## Example Run on Gold-Standard Data
![Alt text](dev_im-3.png)

## Test Run DeepDecode

A proof of concept test bed of DeepDecode is included as a jupyter notebook within this directory. To implement DeepDecode or to request a copy, please reach out to the author directly.

## Author

[@Harshaan Sekhon](https://www.linkedin.com/in/shaan-sekhon-1a217b154/)

## Version History

* 1.0.0: Initial Release

## License

Copyright (C) California Institute of Technology - All Rights Reserved

See LISCENCE.md for more details.
