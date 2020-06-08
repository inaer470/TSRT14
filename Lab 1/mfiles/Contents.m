% Signals and Systems Lab for Matlab
%
% ##
%
%UNDERLYING PRINCIPLES
% The general ideas are summarized below:
%    m=model(description) sets up a model
%    y=simulate(model) simulates a model
%    m=estimate(model,y) estimates a model from data
%    m=rand(model) generates a random object of a specific model
% There are many model classes (nl,lss,ltf,arx,ar,fir,rarx,jarx,ft,covfun,spec,etc)
% but only one signal class (sig).
% Models can be defined by their constructor, by conversion from another model
% or from estimation.
%
%1. Non-linear filtering and estimation
%  nl         Main model object
%  sensormod  Sensor model object, special case of nl
%  sig        Signal object
%  Estimation methods for sensormod
%    ls,wls,ml,crlb
%  Nonlinear filtering methods for nl
%    ekf, ukf, nltf, pf, pmf, crlb
%
%2. CLASSES FOR SIGNAL PROCESSING
%  sig    The signal object
%  ft     The Fourier transform object
%  spec   The spectrum object
%  covfun The covariance function object
%
%  Related m-files: dbsig, getsignal, filtfilt, ncfilter, getfilter
%
%3. CLASSES FOR SYSTEMS
%  lti   Linear Time-Invariant models, parent of ltf and lss
%  ltf   Transfer Functions: continuous or discrete time, SISO or MIMO,
%        certain or uncertain
%  lss   State Space models: deterministic (as TF) or stochastic,
%        continuous or discrete time, SISO or MIMO, certain or uncertain
%  nl    Non-linear time-varying systems, MIMO, certain or uncertain
%  freq  Frequency domain response of LTI models
%
%  Related m-files: exlti, exnl, getfilter
%
%4. CLASSES FOR MODEL ESTIMATION
%  arx   ARX models for discrete time stochastic input-output systems
%  ar    AR models for time-series analysis
%  fir   Finite Impulse Response models
%
%  Related m-files: rand, simulate, estimate
%
%5. CLASSES FOR ADAPTIVE FILTERING AND NON-STATIONARY SIGNAL ANALYSIS
%
%  rarx  Recursive ARX models for time-varying stochastic systems
%  jarx  Jump ARX models for time-varying stochastic systems
%  tfd   Time-Frequency Description using the short-time Fourier Transform
%
%  Related m-files: rand, simulate, estimate
%
%6. CLASSES FOR STATISTICS
%  pdfclass  Parent of all stochastic distribution families:
%  empdist, ndist, expdist, chi2dist, udist, fdist, tdist, betadist, gammadist
%
%  Related m-files: estimate, rand
