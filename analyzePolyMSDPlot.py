#PM env deepmd3

import pandas as pd
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter1d
# from scipy.optimize import curve_fit
from sklearn.linear_model import BayesianRidge
import numpy as np

basicPlot = False
curveFitPlot = False
kernelResample = True

allMSDs = []
# for i in range(3):
for i in [1]:
    df = pd.read_csv(f'{i}MSD.csv')
    allMSDs.append(df)

if kernelResample:

    from sklearn.kernel_ridge import KernelRidge
    from sklearn.model_selection import GridSearchCV
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.gaussian_process.kernels import RBF

    endTime = 1000. # ps
    interval = 0.1 # ps
    lastIndex = int(endTime/interval)

    x_train = np.array(df['dt'][:lastIndex])[::100] # limit to 1ns and then downsample by 100
    y_train = np.array(df['msd_c'][:lastIndex])[::100]
    x_test = np.linspace(0.0, endTime, 100)

    # method from TI paper - credit Sundaraman/Shah
    #Tune model hyperparameters on entire data:
    np.random.seed(0)
    features = x_train[:, None]
    target = y_train
    print('start search')
    # param_grid = {"kernelridge__alpha": np.logspace(-10, 0, 11),
    #                         "kernelridge__kernel": [RBF(l) for l in np.logspace(-3, 3, 10)]}
    param_grid = {"kernelridge__alpha": np.logspace(-10, 0, 5),
                            "kernelridge__kernel": [RBF(l) for l in np.logspace(-3, 3, 5)]}
    grid = GridSearchCV(make_pipeline(StandardScaler(), KernelRidge()),
                                            param_grid=param_grid, verbose=1)
    
    grid.fit(features, target)
    model = grid.best_estimator_
    print('Best parameters:', grid.best_params_)

    #Run model on several sub-samplings of data:
    nRuns = 10 # 1000
    x_test = np.linspace(x_train.min(), x_train.max(), 100)
    predictions = []
    for iRun in range(nRuns):
            print(iRun)
            #Randomly resample data:
            sel = np.sort(np.random.choice(len(x_train), size=len(x_train), replace=True))
            featuresSel = x_train[sel, None]
            targetSel = y_train[sel]
            #Store predictions from subset model:
            model.fit(featuresSel, targetSel)
            predictions.append(model.predict(x_test[:, None]))
    predictions = np.array(predictions)

    #Plot statistics of predictions
    predict_mean = predictions.mean(axis=0)
    predict_err = predictions.std(axis=0)*1.96  # for 95% confidence interval
    predict_lo = predict_mean - predict_err
    predict_hi = predict_mean + predict_err
    plt.fill_between(x_test, predict_lo, predict_hi, color=(1, 0.7, 0.7))
    plt.plot(x_test, predict_mean, 'r')


    plt.tight_layout()
    plt.savefig('fitplotKernelResample.pdf')

if curveFitPlot:
    
    # based on :
    # https://scikit-learn.org/stable/auto_examples/linear_model/plot_bayesian_ridge_curvefit.html
    endTime = 1000. # ps
    interval = 0.1 # ps
    lastIndex = int(endTime/interval)
    
    n_order = 3 #5 #3
    x_train = df['dt'][:lastIndex]
    y_train = df['msd_c'][:lastIndex]
    x_test = np.linspace(0.0, endTime, lastIndex)
    
    X_train = np.vander(x_train, n_order + 1, increasing=True)
    X_test = np.vander(x_test, n_order + 1, increasing=True)
    reg = BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True, max_iter=1000)
    # reg = BayesianRidge(tol=1e-6, compute_score=True)


    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for i, ax in enumerate(axes):
        # Bayesian ridge regression with different initial value pairs
        if i == 0:
            init = [1 / np.var(y_train), 1.0]  # Default values
        elif i == 1:
            init = [1 / np.var(y_train), 1e-3] #[1.0, 1e-3]
            # print(np.var(y_train))
            reg.set_params(alpha_init=init[0], lambda_init=init[1])
        reg.fit(X_train, y_train)
        ymean, ystdRidge = reg.predict(X_test, return_std=True)
        ystd = ystdRidge
        # ax.plot(x_test, func(x_test), color="blue", label="sin($2\\pi x$)")
        ax.scatter(x_train, y_train, s=50, alpha=0.5, label="observation")
        ax.plot(x_test, ymean, color="red", label="predict mean")

        # ystd = np.zeros(lastIndex)
        # for i, yval in enumerate(y_train)

        ax.fill_between(
            x_test, ymean - ystd, ymean + ystd, color="pink", alpha=0.5, label="predict std"
        )
        # ax.set_ylim(-1.3, 1.3)
        ax.legend()
        title = "$\\alpha$_init$={:.2f},\\ \\lambda$_init$={}$".format(init[0], init[1])
        if i == 0:
            title += " (Default)"
        ax.set_title(title, fontsize=12)
        text = "$\\alpha={:.1f}$\n$\\lambda={:.3f}$\n$L={:.1f}$".format(
            reg.alpha_, reg.lambda_, reg.scores_[-1]
        )
        ax.text(0.05, -1.0, text, fontsize=12)

    plt.tight_layout()
    plt.savefig('fitplotBaysRidge.pdf')

if basicPlot:

    # plot all MSDs
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    plt.figure(figsize=(6,4))


    for i in range(3):
        # print(allMSDs[i])
        # plt.loglog(allMSDs[i]['dt'],allMSDs[i]['msd_c'])
        plt.plot(allMSDs[i]['dt'],gaussian_filter1d(allMSDs[i]['msd_c'],3))

    plt.legend(['Pt/O Surface','Bulk','Pt Surface'])  #, loc=2, prop={"size": 20}
    # plt.axis('square')
    # plt.ylim(1,2e3)
    # plt.xlim(1,2e3)
    plt.ylabel("MSD ($\\AA^2$)")
    plt.xlabel("t (ps)")
    # ax = plt.gca()
    # ax.set_aspect('equal', adjustable='box')
    plt.savefig('allMSDs.pdf')
    plt.close()
