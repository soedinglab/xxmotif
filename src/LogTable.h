#ifndef LOG_TABLE_H_
#define LOG_TABLE_H_

class LogTable{
public:
	LogTable();
	~LogTable();

	void setMotifNbCorrection();

	static double      *LOG_i;				/* precalculated logs */
	static double		*LOG_sum_i;			/* precalculated sum of log(1) + .. + log(i) */
	static double 		LOG1_1000;	        /* precalculated log(1.0/1000) */
	static double 		LOG_Neff_pwm;		/* precalculatd log(G->neff_pwm) */
	static double 		LOG_10;				/* precalculated log(10) */
	static double 		LOG_2;				/* precalculated log(2) */
	static double 		LOG_LOG_Seq;		/* precalculated log(log(nSeq(S))) */

	static double*		motifNbCorrection;
};




#endif /* LOG_TABLE_H_ */
