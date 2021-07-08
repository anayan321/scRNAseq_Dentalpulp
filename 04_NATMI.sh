conda create --name NATMI python=3.7.6
conda activate NATMI
sudo apt-get install -y graphviz-dev
pip install pandas=1.0.3 XlsxWriter==1.2.8 xlrd==1.2.0 seaborn==0.10.1 igraph==0.8.3 NetworkX==2.4 PyGraphviz==1.5 bokeh==2.0.2 holoviews==1.13.2
conda install graphviz
git clone https://github.com/asrhou/NATMI.git
cd NATMI

python ExtractEdges.py --species human --emFile all.deep.sc.em.txt --annFile all.deep.sc.ann.txt --interDB lrc2p --coreNum 4 --out ./new_outputs/all.deep.lrc2p.sc
python ExtractEdges.py --species human --emFile all.sound.sc.em.txt --annFile all.sound.sc.ann.txt --interDB lrc2p --coreNum 4 --out ./new_outputs/all.sound.lrc2p.sc

python DiffEdges.py --refFolder new_outputs/all.sound.lrc2p.sc --targetFolder new_outputs/all.deep.lrc2p.sc --interDB lrc2p --out ./new_outputs/DiffEdges_all.sound.deep.lrc2p

python VisInteractions.py --sourceFolder new_outputs/all.sound.lrc2p.sc --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat png --drawNetwork y --plotWidth 8 --plotHeight 8 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6

python VisInteractions.py --sourceFolder new_outputs/all.deep.lrc2p.sc --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat png --drawNetwork y --plotWidth 8 --plotHeight 8 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6

python VisInteractions.py --sourceFolder new_outputs/DiffEdges_all.sound.deep.lrc2p --interDB lrc2p --weightType mean --detectionThreshold 0.2 --plotFormat png --drawNetwork y --plotWidth 10 --plotHeight 10 --layout circle --fontSize 15 --edgeWidth 6 --maxClusterSize 0 --clusterDistance 0.6
