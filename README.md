# Noisy: Predict Natural Noise
Simple scripts for predicting / measuring noise in radio data  and images
(C) Paolo Siera, 2017

```
Usage:
noisy_predictrms.py <file1.ms> [file2.ms [file3.ms ... [fileN.ms] ... ]] [-tsyseff <Tsys/eff (K) OR file>] 
       [-diam <antenna diameter (m)>] [-field <field name>] [-plot <plot name with extension>]

 If you give Tsys/eff as a file it should have two columns: frequency (Hz) and Tsys/eff (K)
```

# Install:
Install through virtual environment:
```
0. Check out from github
1. virtualenv noisy_env
2. source noisy_env/bin/activate
3. pip install [checkout_dir]
4. (noisy_env) noisy_predictrms.py --help
```
