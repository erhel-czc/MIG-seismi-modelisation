import pickle

with open('results.pkl', 'rb') as f:
    result = pickle.load(f)
f.close()

result.slip_rate_evolution()
result.phase_portrait()