function y=dbsignal(name,help)
%DBSIGNAL is the STPK signal database
%   y=dbsignal(name,help)
%   if any argument help is provided, then help text
%   and some introductory plots are shown
%
%   Name      Description
%   bach      A piece of music performed by a cellular phone
%   carpath   Car position obtained by dead-reckoning of wheel velocities
%   current   Current in an overloaded transformator
%   eeg_human The EEG signal y shows the brain activity of a human test subject
%   eeg_rat   The EEG signal y shows the brain activity of a rat.
%   ekg       An EKG signal showing human heart beats.
%   equake    Earthquake data where each of the 14 columns shows one time series.
%   ess       Human speech signal of 's' sound
%   fricest   Data z for a linear regression model used for friction estimation.
%   fuel      Data y=z from measurements of instantanous fuel consumption.
%   genera    The number of genera on earth during 560 million years
%   highway   Measurements of car positions from a helicopter hovering over a highway.
%   pcg       An PCG signal showing human heart beats.
%   planepath Measurements y=p of aircraft position.

% Copyright Fredrik Gustafsson
%$ Revision: 28-Oct-2019 $

y=load([name]);
if nargin>1
    feval(name);
end
