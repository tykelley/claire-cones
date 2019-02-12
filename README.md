This is a repository based on a student in my research group making a comment that the code being used took a few days to run. Based on what was trying to be accomplished, that sounded like an unreasonable amount of time and offered to take a look. We went through the code together in two iterations and I hope to capture and share that interaction because I see these tendencies often.

# Astronomy Background

The code being presented here deals with a quantity in astronomy called the correction factor (which is just a way of saying we know we don't have the right answer, but we think we can get close). This correction factor is correcting our counts of small galaxies in the Milky Way's neighborhood. Smaller galaxies are less bright making them harder to see on the night sky. What this means for observations is that galaxies that are very close tend to be seen easier than those farther away. We can account for this [observational bias](https://en.wikipedia.org/wiki/Malmquist_bias) by looking at numerical simulations of the Milky Way and count the number of galaxies of a given mass for a bunch of different distances from the Milky Way. Comparing the simulations to what we observe results in a plot similar to the one below: 

![Correction Factor](img/corr_factors.png) 

The correction factor is what astronomers use to inform their predictions for the number of galaxies around the Milky Way and, by extension, other galaxies like ours.

# Original Version
