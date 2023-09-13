# Steps for making data vectors

## Baseline:
1- Modify `shear.ini`: under [DEFAULT], change	BASE_DIR	to point to wherever cosmosis-standard-library is (eg. for Lucas this is `/home/secco/cosmosis-standard-library`)

2- Run `cosmosis shear.ini` - this will	create a new file called simulated_shear.fits

3- Modify `shear.ini` again: under [DEFAULT], change existing entries to: 

`2PT_FILE = simulated_shear.fits`

`MODE = 2pt_like`

4- Run `cosmosis shear.ini`, check that	the likelihood is extremely close to zero (eg ~1e-28)

## Contaminated:

### Case 1:
...

### Case 2:
...

