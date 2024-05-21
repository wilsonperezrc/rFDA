<h1>Electronic resources to the article 'Robust Outlier Detection for Functional Data'</h1>

Complementary material to outlier detection in functional data analysis. This approach use the robust Minimum Covariance Determinant estimator to compute the Mahalanobis distance applied to functional principal component scores. 

<h2>R packages:</h2>

Program: R, version 4.0.3

Used R packages:

'fda.usc', version: 6.0-84
'fdaoutlier', version: 1.9.4

<h2>Data sets:</h2>

 We apply the methods analyzed in the study to two real datasets for outlier detection of functional curves.

 <h3>El Niño Sea Surface Temperature data</h3>

 The dataset corresponds to the monthly average of El Niño Sea Surface Temperature (SST), representing Pacific Ocean temperatures in the South American coastal region (coordinates $90^{\circ}$-$80^{\circ}$ West, $0^{\circ}$-$10^{\circ}$ South). To replicate the \textit{rainbow} example \citep{Hyndman_2010}, the dataset to be used is obtained from the \textit{rainbow\_1.4.zip} library, available at \href{http://}{https://cran.r-project.org/src/contrib/Archive/rainbow/}.

  <h3>SCImago Journal Rank data</h3>

A second example, we utilized data from the SCImago Journal Rank website, available at \href{http://}{https://www.scimagojr.com}. This platform annually calculates various metrics to evaluate and rank scientific journals in Scopus. The SCImago system organizes journals into 27 knowledge areas and each thematic area is grouped into specific categories that share common themes, simplifying the search and analysis of scientific publications by thematic field. This study examines the SJR metric values from 1999 to 2022 for the Statistics and Probability category, which is part of the \textit{Mathematics} knowledge area. We focus on journals with SJR values for each year of the time interval. A total of 95 journals satisfy these criteria and are included in this functional study.
