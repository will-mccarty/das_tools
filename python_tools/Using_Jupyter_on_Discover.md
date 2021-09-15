These are my quick notes/instructions for using jupyter notebooks on discover:

* Create a ssh tunnel

First, pick a 5-digit port number.  For this example, I will use 32123, but you should pick one unique to you.  What we need to do is open a port from your machine to discover.  However, you cannot do this in the forward direction.  What we do is open a 'reverse' port forward.  

Mac:  You need to create a 'reverse' port tunnel from discover to your laptop:
```
discover>   ssh -N -f -R 32123:localhost:32123 userid@mylaptop.gsfc.nasa.gov
```
Windows: You cannot ssh directly into a Windows laptop.  Thus, we'll reverse open the port from discover to thunder.gsfc.nasa.gov, and then we will open the port from your laptop to thunder.  This creates a two-hop port forward.
```
discover>   ssh -N -f -R 32123:localhost:32123 userid@thunder.gsfc.nasa.gov
```
Then in putty, add the port forward to 
![putty_screenshot.png](putty_screenshot.png)


