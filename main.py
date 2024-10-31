import numpy as np

from scipy.io import wavfile
import scipy.linalg as la
from itertools import product
import matplotlib.pyplot as plt
import solve_f
import timeit
import tqdm
import ray
from ray.experimental import tqdm_ray
from multiprocessing import cpu_count

np.show_config()
exit(0)

NUM_WORKERS = cpu_count()//2
# ray init can take 3 seconds or more to load
ray.shutdown()

context = ray.init(num_cpus=NUM_WORKERS, include_dashboard=False)
remote_tqdm = ray.remote(tqdm_ray.tqdm)

def time_ifo_N(nFir):
    nbMic = 18
    nbHP = 20
    beta = 0.5
    regu = 1e-0

    Fs, data = wavfile.read('ri_comput.wav')

    # print(data.shape)
    # nFir = data.shape[0] // 10

    HB = np.zeros((nbMic // 2 * (2 * nFir - 1), nbHP * nFir), dtype=np.float32)
    HD = np.zeros((nbMic // 2 * (2 * nFir - 1), nbHP * nFir), dtype=np.float32)
    HB_d = np.zeros((nFir, nbMic), dtype=np.float32)
    # print(HB.shape)
    # exit(0)
    for i, x in enumerate(product(np.arange(nbMic // 2), np.arange(nbHP))):
        HB[x[0] * (2 * nFir - 1):(x[0] + 1) * (2 * nFir - 1), x[1] * nFir:(x[1] + 1) * nFir] = la.toeplitz(
            np.concatenate((data[:nFir, i], np.zeros(nFir - 1))), np.concatenate(([data[0, i]], np.zeros(nFir - 1))))
        HD[x[0] * (2 * nFir - 1):(x[0] + 1) * (2 * nFir - 1), x[1] * nFir:(x[1] + 1) * nFir] = la.toeplitz(
            np.concatenate((data[:nFir, i + 180], np.zeros(nFir - 1))),
            np.concatenate(([data[0, i + 180]], np.zeros(nFir - 1))))
        if x[1] == 4:
            HB_d[:, x[0]] = data[:nFir, i]

    d = np.block(
        [np.concatenate((2e-5 * 10 ** (94 / 20) * HB_d[:, no_Mic], np.zeros(nFir - 1))) for no_Mic in
         range(nbMic // 2)])

    mat_inv = (beta * HB.T @ HB + (1 - beta) * HD.T @ HD + regu * np.eye(HB.shape[1]))
    # print(mat_inv.shape)
    vec = beta * HB.T @ d

    timing = np.zeros(3)
    start = timeit.default_timer()
    np.linalg.inv(mat_inv) @ vec
    timing[0] = (timeit.default_timer() - start)
    start = timeit.default_timer()
    solve_f.solve_f.solvegd(mat_inv.shape[0], beta, HB, d, mat_inv, HB.shape[0], HB.shape[1])
    timing[1] = (timeit.default_timer() - start)
    start = timeit.default_timer()
    solve_f.solve_f.solvegdc(mat_inv.shape[0], beta, HB, d, mat_inv, HB.shape[0], HB.shape[1])
    timing[2] = (timeit.default_timer() - start)

    return timing

# @ray.remote
def time_ifo_N_stat(nFir, N):
    timing = np.zeros((3, N))
    for i in range(N):
        timing[:, i] = time_ifo_N(nFir)

    # bar.update.remote(1)
    return np.concatenate([np.mean(timing, axis=1), np.std(timing, axis=1)])


if __name__ == '__main__':

    nFir_test = np.arange(10, 4000, 20)
    # print(len(nFir_test))

    # bar = remote_tqdm.remote(total=len(nFir_test))
    # stat = np.array(ray.get([time_ifo_N_stat.remote(n, 10, bar) for n in nFir_test]))
    # # bar.close()

    stat = np.zeros((len(nFir_test), 6))
    for i, n in tqdm.tqdm(enumerate(nFir_test), total=len(nFir_test)):
        stat[i, :] = time_ifo_N_stat(n, 5)
        print(time_ifo_N(n))

    plt.figure()
    plt.plot(nFir_test, stat[:, 0], label='Direct inversion')
    plt.fill_between(nFir_test, stat[:, 0] - stat[:, 3], stat[:, 0] + stat[:, 3], alpha=0.2, label='Bande d\'erreur (Variance)')
    plt.plot(nFir_test, stat[:, 1], label='Gradient descent')
    plt.fill_between(nFir_test, stat[:, 1] - stat[:, 4], stat[:, 1] + stat[:, 4], alpha=0.2, label='Bande d\'erreur (Variance)')
    plt.plot(nFir_test, stat[:, 2], label='Conjugate gradient descent')
    plt.fill_between(nFir_test, stat[:, 2] - stat[:, 5], stat[:, 2] + stat[:, 5], alpha=0.2, label='Bande d\'erreur (Variance)')
    plt.xlabel('nFir')
    plt.ylabel('Execution time (s)')
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()

    # start = timeit.default_timer()
    # solve_f.solve_f.solvegd(mat_inv.shape[0], beta, HB, d, mat_inv, HB.shape[0], HB.shape[1])
    # print("PMGd: time = {:.2f} ms".format((timeit.default_timer() - start) * 1000))
    # start = timeit.default_timer()
    # solve_f.solve_f.solvegdc(mat_inv.shape[0], beta, HB, d, mat_inv, HB.shape[0], HB.shape[1])
    # print("PMGdc: time = {:.2f} ms".format((timeit.default_timer() - start) * 1000))
    # start = timeit.default_timer()
    # np.linalg.inv(mat_inv) @ vec
    # print("PM: time = {:.2f} ms".format((timeit.default_timer() - start) * 1000))

    # plt.figure()
    # plt.plot(h)
    # plt.tight_layout()
    # plt.show()

