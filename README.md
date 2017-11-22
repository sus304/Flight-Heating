# Flight-Heating
飛翔中の空力加熱による表面温度を計算するスクリプト。

非常な簡易的なモデルで組み立てられており、飛翔体前面の厚みと物性値、高度履歴と速度履歴のみで計算される。

![hayabusa_qdot](https://user-images.githubusercontent.com/8069773/33128448-cc676a92-cfcf-11e7-9a00-b2bf8a5d6718.png)


対流熱流束はDetra-Kemp-Riddellの方法で計算。
再放射熱流束はTauberの方法で計算。
表面温度は２つの熱流束の和と即座に等しくなるという仮定。

与えられたアブレーション温度に対して表面温度の差が気化熱で吸収できるようにアブレーション厚さを決定する。

## 参考
* 宇宙飛行体の熱気体力学
* Heat Transfer to Satellite Vehicles Re-entering the Atomosphere
* 超軌道速度飛行体の輻射加熱環境に関する研究


# Flight-Heating
This script calculates the surface temperature by aerodynamic heating during flight.

Calculate the object surface tempreture by aerodynamic heating at the time of reentry into the earth, ablation progression during ablation cooling.

Surface Temperature Model: The sum of the heat flux given from the flow and the re-radiation heat is equal to the surface temperature.

Ablation Model: Determine the ablation thickness so that the difference in surface temperature with respect to the set ablation temperature can be eliminaterd by vaporization heat.

## Reference
* 宇宙飛行体の熱気体力学
* Heat Transfer to Satellite Vehicles Re-entering the Atomosphere
* 超軌道速度飛行体の輻射加熱環境に関する研究
