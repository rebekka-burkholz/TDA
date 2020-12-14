# TDA: Tree Distribution Approximation for the Independent Cascade Model (AAAI'21) [[arXiv]](https://arxiv.org/abs/1909.05416)

This repository shares the accompanying code for the publication  "Cascade Size Distributions: Why they matter and How to Compute Them Efficiently" (AAAI 2021) by Rebekka Burkholz and John Quackenbush.

<object data="https://github.com/rebekka-burkholz/TDA/blob/main/fig/Overview.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/rebekka-burkholz/TDA/blob/main/fig/Overview.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>

Efficient algorithms that compute the cascade size distribution and activation probabilities conditional on the final cascade size are available for the Independent Cascade Model (ICM) and implemented in C++, R, and Python.

While a native R implementation is included for readability and easy extensions, we recommend the use of the C++ code, which can be called from R or Python and scales to large networks.
Small tutorials with examples are available in the respective folders.

In case you would like to cite this work, we would like to make it easy for you:
```
@inproceedings{tda,
  title={Cascade Size Distributions: Why they matter and How to Compute Them Efficiently},
  author={Burkholz, Rebekka and Quackenbush, John},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  year={2021},
  organization={AAAI'2021}
}
```
