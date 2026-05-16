import numpy as np


def im2col(I, FH, FW, S=1, P=0):
    N, C, H, W = I.shape
    OH = (H + 2*P - FH) // S + 1
    OW = (W + 2*P - FW) // S + 1
    I = np.pad(I, [(0, 0), (0, 0), (P, P), (P, P)], 'constant')
    C = np.zeros((N, C, FH, FW, OH, OW))
    for y in range(FH):
        ymax = y + S*OH
        for x in range(FW):
            xmax = x + S*OW
            C[:, :, y, x, :, :] = I[:, :, y:ymax:S, x:xmax:S]
    return C.transpose(0, 4, 5, 1, 2, 3).reshape(N*OH*OW, -1)


def col2im(col, input_shape, FH, FW,  S=1, P=0):
    N, C, H, W = input_shape
    OH = (H + 2*P - FH) // S + 1
    OW = (W + 2*P - FW) // S + 1
    col = col.reshape(N, OH, OW, C,  FH, FW).transpose(0, 3, 4, 5, 1, 2)
    I = np.zeros((N, C, H + 2*P + S - 1, W + 2*P + S - 1))
    for y in range(FH):
        ymax = y + S*OH
        for x in range(FW):
            xmax = x + S*OW
            I[:, :, y:ymax:S, x:xmax:S] += col[:, :, y, x, :, :]
    return I[:, :, P:P+H, P:P+W]


def im2colK(I, FH, FW, Sh=1, Sw=1, Ph=0, Pw=0):
    #print(I.shape)
    N, C, H, W = I.shape
    OH = (H + 2*Ph - FH) // Sh + 1
    OW = (W + 2*Pw - FW) // Sw + 1
    #print("OH, OW at im2clK", OH, OW)
    I = np.pad(I, [(0, 0), (0, 0), (Ph, Ph), (Pw, Pw)], 'constant')
    #print(I.shape)
    C = np.zeros((N, C, FH, FW, OH, OW), dtype='float32')
    for y in range(FH):
        ymax = y + Sh*OH
        for x in range(FW):
            xmax = x + Sw*OW
            #print(x,y,I[:, :, y:ymax:Sh, x:xmax:Sw].shape)
            C[:, :, y, x, :, :] = I[:, :, y:ymax:Sh, x:xmax:Sw]
    return C.transpose(0, 4, 5, 1, 2, 3).reshape(N*OH*OW, -1)


def col2imK(col, input_shape, FH, FW, Sh=1, Sw=1, Ph=0, Pw=0):
    N, C, H, W = input_shape
    OH = (H + 2*Ph - FH) // Sh + 1
    OW = (W + 2*Pw - FW) // Sw + 1
    col = col.reshape(N, OH, OW, C,  FH, FW).transpose(0, 3, 4, 5, 1, 2)
    I = np.zeros((N, C, H + 2*Ph + Sh - 1, W + 2*Pw + Sw - 1), dtype='float32')
    for y in range(FH):
        ymax = y + Sh*OH
        for x in range(FW):
            xmax = x + Sw*OW
            I[:, :, y:ymax:Sh, x:xmax:Sw] += col[:, :, y, x, :, :]
    return I[:, :, Ph:Ph+H, Pw:Pw+W]
