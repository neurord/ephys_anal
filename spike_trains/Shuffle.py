from __future__ import print_function as _, division as _
import numpy as np

#different methods of shuffling
def shuffle_everything(T):
    dT = np.diff(T, axis=1)
    dT = np.hstack((T.T[:1].T, dT))
    dT = np.concatenate([dT[None,...] for _ in range(shuffles)])

    assert dT.shape[1] == netsize
    assert dT.shape[0] == shuffles

    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        for s in range(shuffles):
            np.random.shuffle(dT[s, i, :length])

    newT = np.cumsum(dT, axis=2)
    print('newT.shape =', newT.shape, 'newT.dtype =', newT.dtype)
    return newT

def shuffle_bursts(T):
    dT = np.diff(T, axis=1)
    dT = np.hstack((T.T[:1].T, dT))
    dT = np.concatenate([dT[None,...] for _ in range(shuffles)])

    assert dT.shape[1] == netsize
    assert dT.shape[0] == shuffles

    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        for s in range(shuffles):
            dt = dT[s, i, :length]
            indices = np.hstack(([0],
                                 np.where(dt > isiBurstCutoff)[0]))
            bursts = np.hstack((np.diff(indices),
                                [dt.size - indices[-1]]))
            k = np.arange(indices.size)
            np.random.shuffle(k)
            indices = indices[k]
            bursts = bursts[k]

            out = np.empty_like(dt)
            offset = 0
            for j in range(indices.size):
                c = bursts[j]
                out[offset:offset + c] = dt[indices[j]:indices[j] + c]
                offset += c
            dT[s, i, :length] = out

    newT = np.cumsum(dT, axis=2)
    print('newT.shape =', newT.shape, 'newT.dtype =', newT.dtype)
    return newT

def shuffle_everything_first_long(T):
    dT = np.diff(T, axis=1)
    dT = np.hstack((T.T[:1].T, dT))
    dT = np.concatenate([dT[None,...] for _ in range(shuffles)])
    #
    assert dT.shape[1] == netsize
    assert dT.shape[0] == shuffles
    #
    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        any_long = (dT[0, i] > isiBurstCutoff).any()
        for s in range(shuffles):
            while True:
                np.random.shuffle(dT[s, i, :length])
                # if dT[s, i, 0] < isiBurstCutoff and T[i, 1] - T[i, 0] > isiBurstCutoff:
                if any_long and dT[s, i, 0] < isiBurstCutoff:
                    print('reshuffling', i)
                    continue
                break
    #
    newT = np.cumsum(dT, axis=2)
    print('newT.shape =', newT.shape, 'newT.dtype =', newT.dtype)
    return newT

def shuffle_everything_keep_initial(T, hedgehog_mode=False):
    dT = np.diff(T, axis=1)
    dT = np.hstack((T.T[:1].T, dT))
    dT = np.concatenate([dT[None,...] for _ in range(shuffles)])
    #
    assert dT.shape[1] == netsize
    assert dT.shape[0] == shuffles
    #
    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        for s in range(shuffles):
            np.random.shuffle(dT[s, i, 1:length])
    #
    if hedgehog_mode:
        for s in range(shuffles):
            which = -np.isnan(dT[s, :, 0])
            initial = dT[s, which, 0]
            np.random.shuffle(initial)
            dT[s, which, 0] = initial
    #
    newT = np.cumsum(dT, axis=2)
    print('newT.shape =', newT.shape, 'newT.dtype =', newT.dtype)
    return newT

def shuffle_everything_keep_initial_gaussian_shift(T):
    dT = np.diff(T, axis=1)
    dT = np.hstack((T.T[:1].T, dT))
    dT = np.concatenate([dT[None,...] for _ in range(shuffles)])
    #
    assert dT.shape[1] == netsize
    assert dT.shape[0] == shuffles
    #
    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        for s in range(shuffles):
            np.random.shuffle(dT[s, i, 1:length])
    #
    for s in range(shuffles):
        which = -np.isnan(dT[s, :, 0])
        dT[s, which, 0] += np.random.normal(0, firstSpikeShift, which.sum())
    #
    newT = np.cumsum(dT, axis=2)
#    newT[newT >= totalTime] -= 2 * totalTime
#    newT[newT <= totalTime] += 2 * totalTime
#    for s in range(newT.shape[0]):
#        for i in range(newT.shape[1]):
#            np.sort(newT[s, i])

    print('newT.shape =', newT.shape, 'newT.dtype =', newT.dtype)
    return newT

def shuffle_everything_twice(T):
    newT = shuffle_everything_keep_initial(T)
    newT += reshuffleOffset
    newT %= totalTime
    newT = np.sort(newT, axis=2)
    #
    dT = np.diff(newT, axis=2)
    dT = np.concatenate((newT[:, :, :1], dT), axis=2)
    for i in range(dT.shape[1]):
        length = nonnancount(dT[0, i])
        for s in range(shuffles):
            np.random.shuffle(dT[s, i, :length])
    newT = np.cumsum(dT, axis=2)
    #
    return newT
