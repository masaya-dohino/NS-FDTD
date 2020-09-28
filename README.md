# NS-FDTDを用いた髪の毛の構造発色シミュレーション
修士課程での研究


## 背景・目的
リアリティの高い毛髪のCG作成のための研究。毛髪は実は微細な構造を有しており、美しい髪のもつ透明感や艶は微細構造による複雑な光の散乱・干渉がもたらす。しかし、従来のCG映像に使用されているレンダリング手法や理論計算では，単純な事象のみを捉えている場合がほとんどである。ここで、FDTD法と呼ばれる電磁場解析手法を用いて、電磁界空間を差分化して時間領域で電場と磁場を求めることで構造色の解析としてしばしば用いられます。この手法を用いて、CGとして適用できるような髪の毛の光学効果を計算することを考えている。NS-FDTDはより高精度化された電磁場解析アルゴリズムである。

## NS-FDTDについて
従来のFDTD法では差分演算子を用いるため、これによる誤差は避けることができない。そこで、波動方程式の解より導き出される高精度差分法(Non-Standard FDTD)が提案されている。
研究では2次元でのNS-FDTDを実装している。詳細について、先輩の論文より抜粋

![image](https://user-images.githubusercontent.com/57475794/94405983-4c3cfc80-01ac-11eb-9a93-c7a97cf68d5e.png)


![image](https://user-images.githubusercontent.com/57475794/94405606-c456f280-01ab-11eb-9081-97dc93ebddd7.png)



## 遠方界変換について
FDTD シミュレーションにおいて必要となる値は、散乱体付近の微小空間ではなく遠方点における電場と磁場の値である。しかし、計算領域を遠方点まで広げることは不可能である。そこで，近傍界の値から遠方界変換を行い、遠方点における仮想的な値を算出する。遠方界は散乱体を囲む閉曲面 S 上の等価電磁流を積分することによって得られることが知られている。このとき、S 上の等価電磁流には入射波と散乱波の二つの寄与があるが、入射波は閉曲面全体では 0 となるため、結果として散乱波のみの値が得られる。

## シミュレーション空間における毛髪モデルの定義
人間の頭部には直径 50～100µm の毛髪が約 10 万本存在するため、毛髪の外観には繊維集合体としての光学特性が反映される。また、毛髪を構成する組織は様々なダメージを受けることにより、その構造が変化する。これらの構造変化は毛髪繊維の光学特性に影響を及ぼし、毛髪外観の変化として認知される 。ここでは、毛髪の基本的な 3 層構造に分かれており、毛髪は表面を覆ううろこ状のキューティクル、メラニン色素を含む繊維状のコルテックス、毛髪中心部に軸として存在するメデュラの 3 つの組織からなる。


![image](https://user-images.githubusercontent.com/57475794/94406819-91156300-01ad-11eb-97e4-2ffa814c8498.png)
（参考：https://www.kao.com/jp/haircare/hair/1-3/）
## シミュレーション結果
色変換へ
