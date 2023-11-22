# fluidsimulations-demystified
Here I refactored Jos Stam's source code of fluid simulations to use matrices.
In my previous program [Fluid Simulations in C#](https://github.com/jensmalmgren/fluid-simulations-in-csharp), I ported Jos Stam's code from his [paper](https://www.dgp.toronto.edu/public_user/stam/reality/Research/pdf/GDC03.pdf). The program is essentially doing its thing like a shell game.

It could be more enjoyable when you try to understand what is going on. Another aspect of the program is that I would like to make the program multithreaded, and that is impossible when matrices are reused. It is easier to make a multithreaded program if the calculations are confined. Ideally, we need a pipelined calculation flow where each calculation results in a new fresh matrix. Also, this makes it easier to understand what is going on.

I found discrepancies between what Jos claimed is happening in the code and what is happening. In diffusion, he is saying an average is calculated over the four adjacent cells, and that is true, but then he moves the calculation one step forward. In this way, he created a moving average. The previous row is also brought into the calculation in the next row. I removed this behavior and made the averages calculate their values without making moving averages.

In addition, I replaced Jos Stam's linear interpolation with how it is done in this [video](https://youtu.be/qsYE1wMEMPA?si=CUMml4_sguLSmjgD). Not because it is more efficient but because it is better to understand.

The boundary checks can be turned on or off. Less computation is better. This program is less likely to go out of hand. The denominator of 1+4k in the diffusion is acceptable, but it can be even lower.

Enjoy!
