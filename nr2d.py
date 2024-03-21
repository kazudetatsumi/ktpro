import numpy as np, matplotlib.pyplot as plt, argparse
from matplotlib.colors import LogNorm


def plot(data, target, pred, savefig='', time_idx=[0.1, 0.25, 0.6], q_idx=[0.05, 0.5, 0.75]):
    # fig = plt.figure(figsize=(8,8))
    fig, ax = plt.subplots(3,3, figsize=(10,10))

    vmax = np.max([np.max(data), np.max(target), np.max(pred)])
    vmin = max(np.min([np.min(data), np.min(target), np.min(pred)]), 0.1)

    ax[0,0].imshow(data, origin='lower', aspect=0.9, norm=LogNorm(), interpolation='none')
    ax[0,0].set_title('Data')
    ax[0,0].set_ylabel('Time / ch')
    ax[0,0].set_xlabel('q / ch')
    ax[0,1].imshow(target, origin='lower', aspect=0.9, norm=LogNorm(vmax=vmax, vmin=vmin), interpolation='none')
    ax[0,1].set_yticks([])
    ax[0,1].set_title('Ground truth')
    ax[0,1].set_xlabel('q / ch')
    ax[0,2].imshow(pred, origin='lower', aspect=0.9, norm=LogNorm(vmax=vmax, vmin=vmin), interpolation='none')
    ax[0,2].set_yticks([])
    ax[0,2].set_title('Predicted')
    ax[0,2].set_xlabel('q / ch')

    N = 3
    # time_idx = [int((i+1)/N*len(data)/2) for i in range(0,N)]
    time_idx = [int(t*len(data)) for t in time_idx]
    # q_idx = np.linspace(len(data)*0.05, len(data)*0.9, N).astype('int')
    q_idx = [int(t*len(data)) for t in q_idx]
    # print(q_idx)
    for i in range(0,N):
        ax[1,i].plot(data[time_idx[i]], 'o', ms=4, label='Data')
        ax[1,i].plot(target[time_idx[i]], label='Ground truth')
        ax[1,i].plot(pred[time_idx[i]], label='Predicted')
        ax[1,i].set_yscale('log')
        # ax[1,i].set_title(f'time: {time_idx[i]}')
        ax[1,i].set_xlabel('q / ch')
        ax[1,i].legend(frameon=False, fontsize='small', handlelength=1.2, title=f'at t = {time_idx[i]}', alignment='left')
        ax[0,0].axhline(time_idx[i], color='r', linestyle='--', lw=1)
        if i>0:
            ax[1,i].set_yticks([])

        ax[2,i].plot(np.array(data[:,q_idx[i]]), 'o', ms=4, label='Data')
        ax[2,i].plot(np.array(target[:,q_idx[i]]), label='Ground truth')
        ax[2,i].plot(np.array(pred[:,q_idx[i]]), label='Predicted')
        ax[2,i].set_yscale('log')
        # ax[2,i].set_title(f'q: {q_idx[i]}')
        ax[2,i].set_xlabel('Time / ch')
        ax[2,i].legend(frameon=False, fontsize='small', handlelength=1.2, title=f'at q = {q_idx[i]}', alignment='left')
        ax[0,0].axvline(q_idx[i], color='r', linestyle='--', lw=1)
        if i>0:
            ax[2,i].set_yticks([])
    ax[1,0].set_ylabel('Ref. Int. / cts')
    ax[2,0].set_ylabel('Ref. Int. / cts')

    for j in range(1,3):
        yrange = np.array([ax[j,i].get_ylim() for i in range(0,N)])
        for i in range(0,N):
            ax[j,i].set_ylim(np.min(yrange[:,0]), np.max(yrange[:,1]))

    plt.tight_layout()
    if savefig!='':
        if not ('.png' in savefig or '.pdf' in savefig):
            savefig += '.pdf'
        plt.savefig(savefig)
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('num', type=int, help='Data number')
    args = parser.parse_args()

    data = str(args.num).zfill(6) + '_data.npy'
    target = str(args.num).zfill(6) + '_target.npy'

    d = np.load(data)
    t = np.load(target)
    plot(d, t, t)

