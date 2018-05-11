# Violin system-identification estimation experiments

The general idea here is to create a nonlinear acoustic model of a violin with an attached voice-coil exciter.

I am starting with a room-response technique of Farina, which was implemented by Mikhail Baranov (@baranovmv). I hope to add a binary stochastic method of determining the impulse response in addition to the swept-sine technique described in Farina's work. 

http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.33.1614&rep=rep1&type=pdf) by Angelo Farina

![alt text](https://github.com/deganii/ViolinResponse/blob/master/violin-synth.jpg "Experimental Setup")

I also want to add higher-order Volterra terms and explore "deep-learning" system identification techniques in the future.