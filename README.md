https://drive.google.com/open?id=1XyaQYFZI9Z-wLtxpW_WED6iMPiquKGx9

This repository contains implementations of the following results:
1. The valuation of European options in the presence of discrete and proportional dividends, assuming that the underlying follows a Black-Scholes process with deterministic volatility between ex-dividend dates. The valuation can be done by solving backward and forward PDEs, and also by an approximate method extended from Gocsei and Sohel for constant volatility.
2. The calibration of local volatility models with a bootstrapping procedure, utilizing the forward PDE method for European option valuation. With this method, local volatility functions can still be obtained even in the presence of dividends.
3. The calibration of Andersen-Buffum convertible bond model and also the valuation of convertible bonds using this model. The Andersen-Buffum procedure is extended to match the whole equity volatility surface, and is not limited to the ATM volatility as described in their paper.
