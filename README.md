# bookWebsiteMatlabCode

## Supplemental material for the book _Numerical Simulation of Optical Wave Propagation with Examples in MATLAB_

This repository contains MATLAB code used for examples and supporting functions on [schmidtwaveopticsbook.org](https://schmidtwaveopticsbook.org). These tools augment the core simulation framework developed in the book.

## ⚠️ Important Note on Dependencies

The scripts in this repository are designed to work alongside the core functions provided on the disc (or digital download) included with the book.

- **Book Code:** Copyrighted by SPIE. You must have these files in your MATLAB path to run the examples in this repository.
- **Supplemental Code:** The new functions and scripts found _here_ are open-source and intended for community use and modification.

## Reference

Jason D. Schmidt, _Numerical Simulation of Optical Wave Propagation with Examples in MATLAB_, SPIE, Bellingham, WA (2010).

Available from [SPIE Digital Library](https://www.spiedigitallibrary.org/ebooks/PM/Numerical-Simulation-of-Optical-Wave-Propagation-with-Examples-in-MATLAB/eISBN-9780819483270/10.1117/3.866274) and [Amazon](https://www.amazon.com/Numerical-Simulation-Propagation-Examples-Monograph/dp/0819483265/).

## Installation & Usage

1. Clone this repository.
2. Ensure the code from the book's disc is added to your MATLAB path:
   ```matlab
   addpath('/path/to/book_code');
   savepath;
   ```
3. Run any example script located in the /examples directory.

## License

The supplemental code in this repository is licensed under the BSD 3-Clause License. See the LICENSE file for full details.

Note: This license does not apply to the original code provided by SPIE with the purchase of the book.

## Development Quality

![MATLAB CI](https://github.com/jdschmidt-opticalscientist/bookWebsiteMatlabCode/actions/workflows/main.yml/badge.svg?branch=main)

We use GitHub Actions to ensure code quality through automated linting and unit testing.
