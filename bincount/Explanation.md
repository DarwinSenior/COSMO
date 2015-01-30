Explaination of the Bincount and current states.
================

As you probably noticed, this directory is about the Cosmology Object Cluster Correlation Bincounting. People before has done really good job on the bin counting and it is probably great to understand the correlation functions and why bin counting is important for correlation. If you want to know what professor Brunner has done before, I will leads you to [this paper][1]. It basically introduce the $O(nlog(n))$ algorithm for bin counting machanism.

Basically put, if we are given a several number of pairs of points (let us denote them $\phi$ and $\theta$), we will first get the distance between each of the points (using haversine angle, and please google it), like $d_{01}, d_{02}, d_{12}$ so on and so forth. And we put each distance to appropriate bins. For example, bin0 has all distance between 0 and $10^{-5}$ (And we assume each bin has exponetial difference growth). It has been done quite beautifully and you can find the source code in the directories.

Now, the code I did was to find a new direction. We are not using 1 dimensional difference as criteria for bin counting. Instead, we seperate $\phi$ and $\theta$, so that we have two dimensional bins. For example bin00 will contain all the point pairs whose $\phi$ direction differences between 0 to $10^{-5}$ and $\theta$ direction differences between 0 to $10^{-5}$

The naïve way of doing this requires all the bin counting, indicates $O(n^2)$. However, it is very parallelable. I originally planed to parallise code using hadoop cluster. There is a work flow including `cat file | mapper.py | reducer.py | combiner.py`. The relavant file should be `mapper.py reducer.py combiner.py binSortMapReduce.py` It basically chunk the nxn grids into smaller square grids for reducer stage to do substraction and bin counting. If you are interested in the algorithm I am very willing to tell you. So, for this direction, you probably need to work out how to do 3 stage map reduce on hadoop cluster. I feel it is very doable, but I was unfortunately spent too much time on the other idea which are harder to parallise, and left no time for this direction. Feel free to contact me if you want to go on. My email dyue2@illinois.edu

The other way which avoids $O(n^2)$ is to avoid all the explicit calcuation of distance but rather sort and settle the position for each bin boundary using binary search, it requires $O(nlognlogn)$, which has a much bigger advantage for data larger than $10^{6}$, but not so easily parallisable. The relavant file is `binCountSort.py`. The algorithm is fairly complecated and feel free if you want the explaination from me. dyue2@illinois.edu

There are other files that I did experiment with and feel free to see. Some of them are of more naïve implementation and are good for sanaty check and worth looking back if the idea goes to the dead end. 

In the end, please pay attention that the 2d bin counting is an idea that bears some assumption that might not hold true for further investigation. We assume the distance of $\phi$ or longtitute is not related to $\theta$. It does not holds true if there is a big difference between $\theta$ and $\phi$ or just each $\theta$ of the two points. So, when doing further investigation please bear in mind that it is not final at all and so please do not spend too much time on optimsation before coming up with a better model first. That is the mistake I made. And I hope whoever want to investigat this problem, my work sheds some light on this perticular subject. At least my work is not completely wasted.



[1]: [http://www.linuxclustersinstitute.org/conferences/archive/2008/PDF/Dolence_98279.pdf].