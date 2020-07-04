Reinforcement learning (RL) can be useful to solve the traffic assignment problem (TAP), as RL treats drivers as autonomous agents whose goal is to learn to select a route in order to achieve an user equilibrium (UE) situation.
However, the RL task can be hard as it involves hundreds or thousands of agents that are learning simultaneously, thus in a non-stationary environment.
Therefore, our work (see references below) aims at finding metrics that characterize the level of difficulty posed. 
This, on its turn, helps the designer of the learning task to select adequate values for the parameters of the model (e.g., learning rate), as well as the exploration-exploitation strategy.

So far we have proposed two directions for this. The first (Stefanello et al. 2016; Oliveira et al. 2017) 
deals mainly with topological characterization of the network, it considers the demand only in what regards the computation of 
the coupling between shortest paths among the various OD pairs.
The code corresponding to the coupling metric can be found in the route_coupling folder:
https://github.com/maslab-ufrgs/route_coupling

The second metric is based on entropy (Redwan and Bazzan 2020). The code can be found in the entropy_metric folder.

References:
Stefanello, F., da Silva, B. C., and Bazzan, A. L. C. (2016). Using topological statistics to bias and accelerate route choice: preliminary findings in synthetic and real-world road
networks. In Proceedings of Ninth International Workshop on Agents in Traffic and Transportation, pages 1–8, New York, USA.
available at: http://ceur-ws.org/Vol-1678/paper11.pdf

Oliveira, Thiago B. F. and Silva, Bruno C. da and Stefanello F. and Zachow, Arthur 
and Bazzan, Ana L. C., 2017. Extending a Coupling Metric for Characterization of 
Traffic Networks: an Application to the Route Choice Problem. 
Proc. of the 11th Workshop-School on Agents, Enviroments, and Applications (WESAAC 2017). 
(see pag. 56 at http://wesaac2017.c3.furg.br/index.php?Itemid=2344&option=bloco_texto&id_site_componente=3596)

C. Redwan and A. L. C. Bazzan. How hard is for agents to learn the user equilibrium? Characterizing traffic networks by means of entropy.
Submitted. 2020.
