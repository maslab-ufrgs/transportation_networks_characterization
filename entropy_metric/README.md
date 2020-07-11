# Entropy Metric
Python implementation of the entropy metric proposed by Camil

## Requirements
 * [Python 3](https://www.python.org/downloads/)
 * [Python Mathematical Expression Evaluator](https://pypi.python.org/pypi/py_expression_eval)
 * [MSA implementation](https://github.com/maslab-ufrgs/MSA)
 * KSP implementation (use KSP3.py in this repository)
 
  ## Networks
 Available at:  [Networks](https://github.com/maslab-ufrgs/transportation_networks)
 
 ## Usage

```sh
python3 entropy_metric.py [OPTIONS]
```

## Results
Entropy value for each OD pair is printed in console at the end of the iterations

## Options

```
arguments:
  -h, --help            Show help message.
  -f FILE               The network file. 
  -k k_ROUTES           Number of k shortest routes.
  -E EPISODES,          Number of episodes. (default: 1000)
```
## Example

Using [Braess_1_4200_10_c1.net](https://github.com/maslab-ufrgs/transportation_networks/blob/master/Braess%20graphs/Braess_1_4200_10_c1.net) as example, the following call gives us the entropy calculated at the end of 1000 episodes regarding the 3 shortest routes.

```sh
python3 entropy_metric.py -f Braess_1_4200_10_c1.net -k 3 -E 1000
```

One run of the above code returned the entropy value 0.664 to Braess_1_4200_10_c1.net.

## Reference

(Redwan and Bazzan 2020)
