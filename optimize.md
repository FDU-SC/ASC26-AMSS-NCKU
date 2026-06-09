OutBdLow2Hi函数中，由于粗层比细层粒度大，导致细层的数据集中在少数粗层中，此时调用transfer_prolong，计算被主要集中在几个粗层的mpi进程上，导致性能较差。

Restrict 粗->细