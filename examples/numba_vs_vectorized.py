import awkward
import vector
import numpy
import numba
import time


@numba.njit
def compute_mass(photons, n_photons):
    n_events = len(photons)
    out = numpy.empty(n_events, numpy.float64)
    for i in range(n_events):
        total = vector.obj(px=0.0, py=0.0, pz=0.0, E=0.0)
        for j in range(n_photons[i]):
            pho = vector.obj(
                pt = photons[i][j].pt,
                eta = photons[i][j].eta,
                phi = photons[i][j].phi,
                mass = photons[i][j].mass
            )
            total = total + pho
        out[i] = total.mass
    return out

# Same function as above, not compiled with numba
def compute_mass_slow(photons, n_photons):
    n_events = len(photons)
    out = numpy.empty(n_events, numpy.float64)
    for i in range(n_events):
        total = vector.obj(px=0.0, py=0.0, pz=0.0, E=0.0)
        for j in range(n_photons[i]):
            pho = vector.obj(
                pt = photons[i][j].pt,
                eta = photons[i][j].eta,
                phi = photons[i][j].phi,
                mass = photons[i][j].mass
            )
            total = total + pho
        out[i] = total.mass
    return out


time_python_loop = []
time_numba_loop = []
time_vectorized = []

n_events = numpy.logspace(1.5, 5.0, num = 20) 
n_events = n_events.astype(numpy.int64)

for n_event in n_events: 
    photons = vector.awk([[dict({x : numpy.random.normal(0, 1) for x in ("px", "py", "pz")}, E = numpy.random.normal(10, 1)) for inner in range(2)] for outer in range(n_event)])

    n_photons = awkward.num(photons)

    # Compute diphoton mass with plain python loop
    start = time.time()
    mass = compute_mass_slow(photons, n_photons)
    elapsed_python_loop = time.time() - start
    print("Time to compute diphoton mass for %d events (plain python loop): %.6f" % (n_event, elapsed_python_loop))

    # Compute diphoton mass with numba-compiled loop
    mass = compute_mass(photons, n_photons) # run once to compile, then check speed of compiled function
    start = time.time()
    mass = compute_mass(photons, n_photons)
    elapsed_numba_loop = time.time() - start
    print("Time to compute diphoton mass for %d events (numba-compiled loop): %.6f" % (n_event, elapsed_numba_loop))

    # Compute diphoton mass with vectorized operations
    start = time.time()
    mass = (photons[:,0] + photons[:,1]).mass
    elapsed_vectorized = time.time() - start
    print("Time to compute diphoton mass for %d events (vectorized operation): %.6f" % (n_event, elapsed_vectorized))

    time_numba_loop.append(1)
    time_python_loop.append(elapsed_python_loop / elapsed_numba_loop)
    time_vectorized.append(elapsed_vectorized / elapsed_numba_loop)

import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(n_events, time_numba_loop, label = "Numba-compiled loop")
ax1.plot(n_events, time_python_loop, label = "Plain python loop")
ax1.plot(n_events, time_vectorized, label = "vectorized operations")

ax1.set_xscale("log")
ax1.set_yscale("log")

plt.xlabel("Number of events")
plt.ylabel("Time (relative to numba-compiled loop)")
plt.legend()

plt.savefig("vector_performance.pdf")
