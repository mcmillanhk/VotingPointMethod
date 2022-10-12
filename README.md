# Voting Point Method

This repository contains a Matlab version of the Voting Point Method code describe in this paper:

McMillan, H.K. and Westerberg, I.K., 2015. Rating curve estimation under epistemic uncertainty. Hydrological Processes, 29(7), pp.1873-1882.

The code version available here was recoded using a simpler structure by David Ocio Moreno, and tested for accuracy against the original version, for this paper:

Ocio, D., Le Vine, N., Westerberg, I., Pappenberger, F. and Buytaert, W., 2017. The role of rating curve uncertainty in real‐time flood forecasting. Water Resources Research, 53(5), pp.4197-4213.

Please cite both papers if using this code.


Note two small corrections to the text in McMillan and Westerberg (2015) as follows:

On page 1877 in the second column and the last paragraph Rg should be the ratio between rated
and gauged discharge for point g.

On page 1878 in the first column and the first paragraph, the logs of h, q, hfit and qfit were taken
for all sections except the uppermost.

# Run the code

To run the code with the example data, use command VPlogistic_3seg('VPMratingdata')
