# -*- coding: utf-8 -*-
"""
Created on Mon May  8 16:18:09 2023

@author: lhh
"""

"""
#安装模块
pip install mrmr_selection
import mrmr

"""

import numpy as np
import pandas as pd
import pickle
import xgboost as xgb
from sklearn.model_selection import KFold, GridSearchCV,cross_val_score,train_test_split, LeaveOneOut
from sklearn.metrics import precision_recall_curve,auc,roc_curve,accuracy_score,confusion_matrix
import matplotlib.pyplot as plt

import os
os.chdir('E:/work/240624_brca/25_shap/genePair_based/result')

seed=123
np.random.seed(seed)

# 准备数据
import pyreadr
featureMatrix = pyreadr.read_r("E:/work/240624_brca/05_sdcn/3000/train/001/data/train_2353_1631_3984.rda")['train3984']
sampleLabel = pd.read_table("E:/work/240624_brca/05_sdcn/3000/train/001/result/3984/result/sdcn/pred2.txt", header = None)
sampleLabel.index = featureMatrix.index
sampleLabel.rename = "cluster"
sampleLabel = sampleLabel.astype(int)
sampleLabel = sampleLabel.iloc[:,0]
# sampleLabel = sampleLabel.replace(to_replace=[0, 1], value=[1, 2])


i = 3984
fpN = featureMatrix.iloc[:, 0:i]
X = fpN
y = sampleLabel


########调整参数
#将弱分类器集成在一起的强分类器
xgb_model=xgb.XGBClassifier(learning_rate=0.1,n_estimators=100,max_depth=3,colsample_bytree=0.8,objective='binary:logistic')

#超参数
parameters={'learning_rate':[0.1,0.01,0.001],'n_estimators':[50,100,200],'max_depth':[3,6,9],'colsample_bytree':[0.6,0.8,1.0]}

#参数寻优
clf=GridSearchCV(xgb_model, parameters, cv=10, scoring='accuracy') #参数调优，寻找最优的模型参数组合
clf.fit(X, y) #训练模型

print("Best Params:",clf.best_params_)
print("Validation Accuracy:",clf.best_score_)

clf.best_params_
clf.best_estimator_

#######交叉验证
# np.random.seed(seed)
xgb_model=clf.best_estimator_

kf = KFold(n_splits=10, shuffle=True, random_state=seed)

#交叉验证
results=cross_val_score(xgb_model, X, y, cv=kf) #评估模型性能
print('Standardize: %.3f (%.3f) MSE' % (results.mean(),results.std()))

print(results.mean())
result1=pd.DataFrame(results)
result_mean=results.mean()
result1.loc[len(result1.index)] = result_mean 
indexName = list(range(1, 11))
indexName.append('mean')
result1.index = indexName
result1.columns = ['cross_val_score']
result1.to_csv('cross_val_score_xgb_10cv.tsv',sep='\t')

xgb_model.fit(X, y)

##########绘制ROC曲线
y_train_prob=xgb_model.predict_proba(X)
fpr,tpr,thresholds=roc_curve(y,y_train_prob[:,1]) #sklearn库计算ROC曲线面积
roc_auc=auc(fpr,tpr) 
y_train_prob = pd.DataFrame(y_train_prob)
y_train_prob.index = X.index
y_train_prob.columns=[1, 2]
y_train_prob.to_csv("train_predict_prob_xgb_10cv.tsv", sep = "\t")

parameters = {'axes.labelsize': 15,
              'axes.titlesize': 15,
              'xtick.labelsize': 13,
              'ytick.labelsize': 13,
              'font.family': 'Times New Roman'
              }
plt.rcParams.update(parameters)


axes=plt.subplots(1,1,figsize=(6,6),dpi=300) #matplotlib库函数，用于创建多个子图
lab='Overall.AUC=%.4f' % (roc_auc)

axes[1].step(fpr,tpr,label=lab,lw=2,color='red')
axes[1].set_title('ROC Curve of Training Set',fontsize=12,fontname='Times New Roman')
axes[1].set_xlabel('False Positive Rate',fontsize=12,fontname='Times New Roman')
axes[1].set_ylabel('True Positive Rate',fontsize=12,fontname='Times New Roman')
axes[1].legend(loc='lower right',prop={'family':'Times New Roman','size':12})
plt.show()

axes[0].savefig('train_roc_xgb_10cv.pdf')

##train数据集上混淆矩阵
import seaborn as sns

y_train_pred=xgb_model.predict(X)
y1 = y
y1 = y.replace(to_replace=[0, 1], value=[1, 2])
y_train_pred[y_train_pred==1]=2
y_train_pred[y_train_pred==0]=1
cm_train=confusion_matrix(y1, y_train_pred)

plt.figure(figsize=(8, 8))
sns.heatmap(cm_train, annot=True,fmt='.20g',cmap='Blues', yticklabels=[1,2], xticklabels=[1,2])#添加fmt，不使用科学计数法显示
plt.title("Confusion Matrix of Training Set")
plt.xlabel('Predicted labels')
plt.ylabel('True labels')
plt.savefig('confusion_matrix_train_xgb_10cv.pdf')
plt.show()

accuracy1 = accuracy_score(y1, y_train_pred)
file1 = open('./train_accuracy_xgb_10cv.txt','w')
file1.write(str(accuracy1))
file1.close()

y_train_pred = pd.DataFrame(y_train_pred)
y_train_pred.index = X.index
y_train_pred.to_csv('train_pred_xgb_10cv.tsv',sep='\t')



######test-ROC

X_test = pyreadr.read_r("E:/work/240624_brca/05_sdcn/3000/test/001/data/test_2353_1631_3984.rda")['test3984']
y_test = pd.read_table("E:/work/240624_brca/05_sdcn/3000/test/001/result/3984/result/sdcn/pred2.txt", header = None)
# y_test.loc[y_test[0]==1]=2
# y_test.loc[y_test[0]==0]=1

y_test.index = X_test.index
y_test.rename = "cluster"

y_test = y_test.astype(int)
y_test = y_test.iloc[:,0]



y_test_prob=xgb_model.predict_proba(X_test)
fpr,tpr,thresholds=roc_curve(y_test,y_test_prob[:,1])
roc_auc=auc(fpr,tpr)
y_test_prob = pd.DataFrame(y_test_prob)
y_test_prob.index = X_test.index
y_test_prob.columns=[1, 2]
y_test_prob.to_csv("test_predict_prob_xgb_10cv.tsv", sep = "\t")

axes=plt.subplots(1,1,figsize=(6,6),dpi=300)
lab='Overall.AUC=%.4f' % (roc_auc)

axes[1].step(fpr,tpr,label=lab,lw=2,color='red')
axes[1].set_title('ROC Curve of Testing Set',fontsize=12,fontname='Times New Roman')
axes[1].set_xlabel('False Positive Rate',fontsize=12,fontname='Times New Roman')
axes[1].set_ylabel('True Positive Rate',fontsize=12,fontname='Times New Roman')
axes[1].legend(loc='lower right',prop={'family':'Times New Roman','size':12})
plt.show()
axes[0].savefig('test_roc_xgb_10cv.pdf')


###################confusion matrix
#可视化地评估监督学习算法的性能（在整个训练集上运用）
import seaborn as sns

y_test_pred=xgb_model.predict(X_test)
y_test1 = y_test
y_test1 = y_test.replace(to_replace=[0, 1], value=[1, 2])
y_test_pred[y_test_pred==1]=2
y_test_pred[y_test_pred==0]=1
cm=confusion_matrix(y_test1, y_test_pred)

plt.figure(figsize=(8, 8))
sns.heatmap(cm, annot=True,fmt='.20g',cmap='Blues', yticklabels=[1,2], xticklabels=[1,2])#添加fmt，不使用科学计数法显示
plt.title("Confusion Matrix of Testing Set")
plt.xlabel('Predicted labels')
plt.ylabel('True labels')
plt.savefig('confusion_matrix_test_xgb_10cv.pdf')
plt.show()

accuracy2 = accuracy_score(y_test1, y_test_pred)
file2 = open('./test_accuracy_xgb_10cv.txt','w')
file2.write(str(accuracy2))
file2.close()

y_test_pred = pd.DataFrame(y_test_pred)
y_test_pred.index = X_test.index
y_test_pred.to_csv('test_pred_xgb_10cv.tsv',sep='\t')




#########XGBoost特征重要性
feature_names=X.columns
feature_importances=xgb_model.feature_importances_
indices = np.argsort(feature_importances)[::-1] #argsort返回的是数组值从小到大的索引值；倒序操作（从大到小的索引值）
for index in indices:
    print("特征 %s 重要度为 %f" %(feature_names[index], feature_importances[index]))

#输出特征重要性
importances=pd.DataFrame(feature_importances)
importances.index=feature_names
importances.columns=['feature_importances']
importances_order=importances.sort_values(by=['feature_importances'],ascending=False) #排序，降序排列
importances_order.to_csv('feature importances.tsv',sep='\t')

# plt.figure(figsize=(35,8))
# plt.title("Feature Importances-" + str(len(feature_importances)), fontsize=30,fontname='Times New Roman')
# plt.bar(range(len(feature_importances)), feature_importances[indices], color='b')
# plt.xticks(range(len(feature_importances)), np.array(feature_names)[indices], rotation=90,color='b',fontsize=16) #获取或设置当前x轴刻度位置和标签
# plt.savefig('feature importances-'+str(len(feature_importances))+'.pdf',bbox_inches='tight')

plt.figure(figsize=(16,8))
plt.title("Feature Importances-20",fontsize=30,fontname='Times New Roman')
plt.bar(range(len(feature_importances))[0:21], feature_importances[indices][0:21], color='#F78179')
plt.xticks(range(len(feature_importances))[0:21], np.array(feature_names)[indices][0:21], rotation=90,color='black',fontsize=20,fontname='Times New Roman')
plt.yticks(fontsize=20,fontname='Times New Roman')
plt.savefig('feature importances-20.pdf',bbox_inches='tight')



###############shap
import shap
shap.initjs() #为了能够输出shap的图像
explainer=shap.Explainer(xgb_model) #解释机器学习模型
shap_values=explainer(X) #为模型创建可解释性对象

shap.plots.force(shap_values)
shap.save_html('shap-feature-all.html',shap.plots.force(shap_values))

shap.plots.force(shap_values[0],matplotlib=True,show=True, text_rotation=5)
shap.plots.force(shap_values[0],matplotlib=True,show=True, text_rotation=5)
shap.plots.force(shap_values[0],matplotlib=True,show=False, text_rotation=5)#这是两个不同分类的特征贡献的两个例子
plt.savefig('shap-feature-single.pdf',bbox_inches='tight')
shap.plots.force(shap_values[1],matplotlib=True,show=True, text_rotation=5)
shap.plots.force(shap_values[1],matplotlib=True,show=True, text_rotation=5)
shap.plots.force(shap_values[1],matplotlib=True,show=False, text_rotation=5)
plt.savefig('shap-feature-single-1.pdf',bbox_inches='tight')


shap.plots.beeswarm(shap_values,max_display=500)#手动保存
#plt.savefig("beeswarm-500.png", dpi = 300)
shap.plots.beeswarm(shap_values,max_display=21)#手动保存
#plt.savefig("beeswarm-20.png", dpi = 300)

shap.plots.bar(shap_values,max_display=500,show=True)#先使用True能够在plot界面画出正确图像，然后show改为False保存
shap.plots.bar(shap_values,max_display=500,show=True)
shap.plots.bar(shap_values,max_display=500,show=False)#show=False可以保存图片
plt.savefig('shap-plot-bar.pdf',bbox_inches='tight')

shap.plots.bar(shap_values,max_display=21,show=True)
shap.plots.bar(shap_values,max_display=21,show=True)
shap.plots.bar(shap_values,max_display=21,show=False)#show=False可以保存图片
plt.savefig('shap-plot-bar-20.pdf',bbox_inches='tight')

shap.plots.waterfall(shap_values[0],max_display=500,show=True)
shap.plots.waterfall(shap_values[0],max_display=500,show=True)
shap.plots.waterfall(shap_values[0],max_display=500,show=False)#两个不同分类的特征贡献的例子
plt.savefig('shap-plot-waterfall-single.pdf',bbox_inches='tight')
shap.plots.waterfall(shap_values[1],max_display=500,show=True)
shap.plots.waterfall(shap_values[1],max_display=500,show=True)
shap.plots.waterfall(shap_values[1],max_display=500,show=False)
plt.savefig('shap-plot-waterfall-single-1.pdf',bbox_inches='tight')

shap.plots.waterfall(shap_values[0],max_display=21,show=True)
shap.plots.waterfall(shap_values[0],max_display=21,show=True)
shap.plots.waterfall(shap_values[0],max_display=21,show=False)#两个不同分类的特征贡献的例子，仅取前20的特征
plt.savefig('shap-plot-waterfall-20-single.pdf',bbox_inches='tight')
shap.plots.waterfall(shap_values[1],max_display=21,show=True)
shap.plots.waterfall(shap_values[1],max_display=21,show=True)
shap.plots.waterfall(shap_values[1],max_display=21,show=False)
plt.savefig('shap-plot-waterfall-20-single-1.pdf',bbox_inches='tight')


## 提取特征并按重要性排序
# 计算每个特征的重要性 (绝对值平均)
shap_importance = pd.DataFrame({
    'feature': X.columns,
    'importance': np.abs(shap_values.values).mean(axis=0)
})

# 按重要性排序
shap_importance = shap_importance.sort_values(by='importance', ascending=False)

# 保存特征到 CSV 文件
shap_importance.to_csv('E:/work/240624_brca/25_shap/genePair_based/data/shap_feature_importance.csv', index=False)
