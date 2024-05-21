<h1>Electronic resources to the article 'Robust Outlier Detection for Functional Data'</h1>

Complementary material to outlier detection in functional data analysis. This approach use the robust Minimum Covariance Determinant estimator to compute the Mahalanobis distance applied to functional principal component scores. 

<h2>R packages:</h2>

Program: R, version 4.0.3

Used R packages:
<UL>
<LI> 'fda.usc', version: 2.0.2
<LI> 'fdaoutlier', version: 0.2.0
<LI> 'raibow', version: 3.7
</UL>


<h2>Data sets:</h2>

 We apply the methods analyzed in the study to two real datasets for outlier detection of functional curves.

 <h3>El Niño Sea Surface Temperature data</h3>

 The dataset corresponds to the monthly average of El Niño Sea Surface Temperature, representing Pacific Ocean temperatures in the South American coastal region (coordinates $90^{\circ} - 80^{\circ}$ West, $0^{\circ} - 10^{\circ}$ South). The dataset to be used is obtained from the <em>rainbow\_1.4.zip</em> library, available at https://cran.r-project.org/src/contrib/Archive/rainbow/.

  <h3>SCImago Journal Rank data</h3>

A second example, we utilized data from the SCImago Journal Rank website, available at https://www.scimagojr.com. This platform annually calculates various metrics to evaluate and rank scientific journals in Scopus. The SCImago system organizes journals into 27 knowledge areas and each thematic area is grouped into specific categories that share common themes. This study examines the SJR metric values from 1999 to 2022 for the <em>Statistics and Probability</em> category, which is part of the <em>Mathematics</em> knowledge area. We focus on journals with SJR values for each year of the time interval. A total of 95 journals satisfy these criteria and are included in this functional study.
