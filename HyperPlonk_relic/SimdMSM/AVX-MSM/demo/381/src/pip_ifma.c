#define _POSIX_C_SOURCE 199309L
#include "pip_ifma.h"
#include"time.h"

double get_elapsed_time(struct timespec start, struct timespec end) {
    return (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
}

void simdPADDaff(uint64_t u1x[][HT_NWORDS], uint64_t u1y[][HT_NWORDS], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], point52_t p3[8])
{
    htpoint a, b, c;
    int i;

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(u1x[7][i], u1x[6][i], u1x[5][i], u1x[4][i], u1x[3][i], u1x[2][i], u1x[1][i], u1x[0][i]);
        a.y[i] = set_vector(u1y[7][i], u1y[6][i], u1y[5][i], u1y[4][i], u1y[3][i], u1y[2][i], u1y[1][i], u1y[0][i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
    }

    gfp_mont_relic2avx_8x1w(a.x, a.x);
    gfp_mont_relic2avx_8x1w(a.y, a.y);
    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.z[i] = VSET1(ht_montR[i]);
        b.z[i] = VSET1(ht_montR[i]);
    }
    PADD_projc_8x1w_z1eqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->x, c.x, 0);
    get_channel_8x1w(p3[1]->x, c.x, 1);
    get_channel_8x1w(p3[2]->x, c.x, 2);
    get_channel_8x1w(p3[3]->x, c.x, 3);
    get_channel_8x1w(p3[4]->x, c.x, 4);
    get_channel_8x1w(p3[5]->x, c.x, 5);
    get_channel_8x1w(p3[6]->x, c.x, 6);
    get_channel_8x1w(p3[7]->x, c.x, 7);

    get_channel_8x1w(p3[0]->y, c.y, 0);
    get_channel_8x1w(p3[1]->y, c.y, 1);
    get_channel_8x1w(p3[2]->y, c.y, 2);
    get_channel_8x1w(p3[3]->y, c.y, 3);
    get_channel_8x1w(p3[4]->y, c.y, 4);
    get_channel_8x1w(p3[5]->y, c.y, 5);
    get_channel_8x1w(p3[6]->y, c.y, 6);
    get_channel_8x1w(p3[7]->y, c.y, 7);

    get_channel_8x1w(p3[0]->z, c.z, 0);
    get_channel_8x1w(p3[1]->z, c.z, 1);
    get_channel_8x1w(p3[2]->z, c.z, 2);
    get_channel_8x1w(p3[3]->z, c.z, 3);
    get_channel_8x1w(p3[4]->z, c.z, 4);
    get_channel_8x1w(p3[5]->z, c.z, 5);
    get_channel_8x1w(p3[6]->z, c.z, 6);
    get_channel_8x1w(p3[7]->z, c.z, 7);
    
}

void simdPADDprj(buf2_st* buf, point52_t p3[8])
{
    htpoint a, b, c;
    int i;

    for(i = 1; i < HT_NWORDS; i++)
    {
        if(buf[i].op1 == NULL)
        {
            buf[i].op1 = buf[0].op1;
            buf[i].op2 = buf[0].op2;
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector((*buf[7].op1)->x[i], (*buf[6].op1)->x[i], (*buf[5].op1)->x[i], (*buf[4].op1)->x[i], (*buf[3].op1)->x[i], (*buf[2].op1)->x[i], (*buf[1].op1)->x[i], (*buf[0].op1)->x[i]);
        a.y[i] = set_vector((*buf[7].op1)->y[i], (*buf[6].op1)->y[i], (*buf[5].op1)->y[i], (*buf[4].op1)->y[i], (*buf[3].op1)->y[i], (*buf[2].op1)->y[i], (*buf[1].op1)->y[i], (*buf[0].op1)->y[i]);
        a.z[i] = set_vector((*buf[7].op1)->z[i], (*buf[6].op1)->z[i], (*buf[5].op1)->z[i], (*buf[4].op1)->z[i], (*buf[3].op1)->z[i], (*buf[2].op1)->z[i], (*buf[1].op1)->z[i], (*buf[0].op1)->z[i]);
        b.x[i] = set_vector((*buf[7].op2)->x[i], (*buf[6].op2)->x[i], (*buf[5].op2)->x[i], (*buf[4].op2)->x[i], (*buf[3].op2)->x[i], (*buf[2].op2)->x[i], (*buf[1].op2)->x[i], (*buf[0].op2)->x[i]);
        b.y[i] = set_vector((*buf[7].op2)->y[i], (*buf[6].op2)->y[i], (*buf[5].op2)->y[i], (*buf[4].op2)->y[i], (*buf[3].op2)->y[i], (*buf[2].op2)->y[i], (*buf[1].op2)->y[i], (*buf[0].op2)->y[i]);
        b.z[i] = set_vector((*buf[7].op2)->z[i], (*buf[6].op2)->z[i], (*buf[5].op2)->z[i], (*buf[4].op2)->z[i], (*buf[3].op2)->z[i], (*buf[2].op2)->z[i], (*buf[1].op2)->z[i], (*buf[0].op2)->z[i]);
    }


    PADD_projc_8x1w_z1neqz2(&c, &a, &b);

    get_channel_8x1w(p3[0]->x, c.x, 0);
    get_channel_8x1w(p3[1]->x, c.x, 1);
    get_channel_8x1w(p3[2]->x, c.x, 2);
    get_channel_8x1w(p3[3]->x, c.x, 3);
    get_channel_8x1w(p3[4]->x, c.x, 4);
    get_channel_8x1w(p3[5]->x, c.x, 5);
    get_channel_8x1w(p3[6]->x, c.x, 6);
    get_channel_8x1w(p3[7]->x, c.x, 7);

    get_channel_8x1w(p3[0]->y, c.y, 0);
    get_channel_8x1w(p3[1]->y, c.y, 1);
    get_channel_8x1w(p3[2]->y, c.y, 2);
    get_channel_8x1w(p3[3]->y, c.y, 3);
    get_channel_8x1w(p3[4]->y, c.y, 4);
    get_channel_8x1w(p3[5]->y, c.y, 5);
    get_channel_8x1w(p3[6]->y, c.y, 6);
    get_channel_8x1w(p3[7]->y, c.y, 7);

    get_channel_8x1w(p3[0]->z, c.z, 0);
    get_channel_8x1w(p3[1]->z, c.z, 1);
    get_channel_8x1w(p3[2]->z, c.z, 2);
    get_channel_8x1w(p3[3]->z, c.z, 3);
    get_channel_8x1w(p3[4]->z, c.z, 4);
    get_channel_8x1w(p3[5]->z, c.z, 5);
    get_channel_8x1w(p3[6]->z, c.z, 6);
    get_channel_8x1w(p3[7]->z, c.z, 7);
    
}

void simdPADDmix(point52 *bucket[AVX_WAY], uint64_t u2x[][HT_NWORDS], uint64_t u2y[][HT_NWORDS], int num)
{
    htpoint a, b, c;
    int i;

    for(i = num; i < HT_NWORDS; i++)
    {
        if(bucket[i] == NULL)
        {
            bucket[i] = bucket[0];
        }
    }

    for (i = 0; i < HT_NWORDS; i++)
    {
        a.x[i] = set_vector(bucket[7]->x[i], bucket[6]->x[i], bucket[5]->x[i], bucket[4]->x[i], bucket[3]->x[i], bucket[2]->x[i], bucket[1]->x[i], bucket[0]->x[i]);
        a.y[i] = set_vector(bucket[7]->y[i], bucket[6]->y[i], bucket[5]->y[i], bucket[4]->y[i], bucket[3]->y[i], bucket[2]->y[i], bucket[1]->y[i], bucket[0]->y[i]);
        a.z[i] = set_vector(bucket[7]->z[i], bucket[6]->z[i], bucket[5]->z[i], bucket[4]->z[i], bucket[3]->z[i], bucket[2]->z[i], bucket[1]->z[i], bucket[0]->z[i]);
        b.x[i] = set_vector(u2x[7][i], u2x[6][i], u2x[5][i], u2x[4][i], u2x[3][i], u2x[2][i], u2x[1][i], u2x[0][i]);
        b.y[i] = set_vector(u2y[7][i], u2y[6][i], u2y[5][i], u2y[4][i], u2y[3][i], u2y[2][i], u2y[1][i], u2y[0][i]);
        b.z[i] = VSET1(ht_montR[i]);
   }

    gfp_mont_relic2avx_8x1w(b.x, b.x);
    gfp_mont_relic2avx_8x1w(b.y, b.y);
    PADD_projc_8x1w_mix(&c, &a, &b);

    for (i = 0; i < num ; i++)
    {
        get_channel_8x1w(bucket[i]->x, c.x, i);
        get_channel_8x1w(bucket[i]->y, c.y, i);
        get_channel_8x1w(bucket[i]->z, c.z, i);
    }
    
}


void avx2_withIFMA(pip_ctx *ctx, int eff_times)
{
    int i;
    buf2_st *item;
    point52_t avx2_out[AVX_WAY];
    
    /* avx2*/
    simdPADDprj(ctx->buf2, avx2_out);

    for (i = 0; i < eff_times; i++)
    {
        item = &(ctx->buf2[i]);
        if(item->flag == 0)
		{
			queue_in(&(ctx->T), avx2_out[i], item->id);
		}
		else
		{
			// ep_copy(ctx->bucket[item->id].B, avx2_out[i]);
            point52_copy(*item->op2, avx2_out[i]);
		}
        
    }

    ctx->cnt2 = 0;
    
}
double total_time_avx2 = 0.0;
/* process of T*/
void state_trans_T_withIFMA(pip_ctx *ctx)
{
	buf2_st *same_id;
	int bucket_id, i;
    struct timespec t_start, t_end;

	while (!queue_is_empty(&(ctx->T))) 
	{
		i = ctx->T.front;
        bucket_id = ctx->T.pbuf[i].id;
		same_id = check_eq_id(ctx, bucket_id);
		if(same_id == NULL)
		{	

            /* write into bucket*/
            ctx->buf2[ctx->cnt2].id = bucket_id;
            ctx->buf2[ctx->cnt2].op1 = &(ctx->T.pbuf[i].point);
            ctx->buf2[ctx->cnt2].op2 = &(ctx->bucket[bucket_id].B);
            ctx->buf2[ctx->cnt2].flag = 1;
            ctx->cnt2++;

		}
		else
		{
			/* write into queue*/
	    	same_id->op2 = &(ctx->T.pbuf[i].point);
			same_id->flag = 0;
		}

		if(ctx->cnt2 == AVX_WAY)
		{
            clock_gettime(CLOCK_MONOTONIC, &t_start);

			avx2_withIFMA(ctx, AVX_WAY);

            clock_gettime(CLOCK_MONOTONIC, &t_end);
            total_time_avx2 += get_elapsed_time(t_start, t_end);
		}

		queue_out(&(ctx->T), 1);
    }
}

int point52_is_infty(const point52_t p) {
	return (fp_is_zero(p->z) == 1);
}

double total_avx_prepare = 0.0; // 阶段1: 格式转换
double total_avx_calc    = 0.0; // 阶段2: SIMD计算
double total_avx_queue   = 0.0; // 阶段3: 结果处理
double total_avx_trans   = 0.0; // 阶段4: 状态转移 (你重点关注的)
void mpi_conv_64to52_8x1w(htfe_t r,const htfe_t a,int rlen,int alen){
    int i = 0, j = 0;
    int shr_pos = 64;
    int shl_pos = 0;

    // 1. 初始化 AVX 常量和变量
    __m512i v_temp = _mm512_setzero_si512(); 
    __m512i v_mask = _mm512_set1_epi64(HT_BMASK);
    __m512i v_current;

    // 2. 预加载第一个输入向量 a[0]
    // 注意：a 是 htfe_t 类型 (即 __m512i 数组)，直接索引即可
    if (j < alen) {
        v_current = a[j];
    } else {
        v_current = _mm512_setzero_si512();
    }

    // 3. 主循环: 按照输出肢的个数迭代 (流式处理)
    // 逻辑完全复刻标量版的 while 循环，但变量全是向量
    while (i < rlen && j < alen)
    {
        // ---------------------------------------------------------
        // A. 并行移位与拼接
        // word = ((temp >> shr_pos) | (a[j] << shl_pos));
        // ---------------------------------------------------------
        
        // 准备移位计数器 (__m128i 格式)
        __m128i k_shr = _mm_cvtsi32_si128(shr_pos);
        __m128i k_shl = _mm_cvtsi32_si128(shl_pos);

        // 执行并行移位
        // 注意：AVX-512 如果移位量 >= 64 会自动置0，这比 C 语言更安全
        __m512i v_right = _mm512_srl_epi64(v_temp, k_shr);
        __m512i v_left  = _mm512_sll_epi64(v_current, k_shl);
        
        // 拼接
        __m512i v_word = _mm512_or_si512(v_right, v_left);
        
        // ---------------------------------------------------------
        // B. 掩码提取并存储
        // r[i] = word & mask;
        // ---------------------------------------------------------
        r[i] = _mm512_and_si512(v_word, v_mask);

        // ---------------------------------------------------------
        // C. 更新状态变量
        // ---------------------------------------------------------
        shr_pos -= 12;
        shl_pos += 12;

        // 判断是否需要消耗下一个输入 Limb
        // 逻辑: if ((shr_pos > 0) && (shl_pos < 64))
        if ((shr_pos > 0) && (shl_pos < 64))
        {
            v_temp = v_current; // 当前输入变为"旧数据"
            j++;
            
            // 加载下一个输入向量
            if (j < alen) {
                v_current = a[j];
            } else {
                v_current = _mm512_setzero_si512();
            }
        }

        // 处理位移量的回绕 (Wrap around)
        if (shr_pos <= 0) shr_pos += 64;
        if (shl_pos >= 64) shl_pos -= 64;
        
        // 特殊情况：如果 temp 恰好被完全消耗完 (例如 64->32 转换时)，重置 temp
        if (shr_pos == 64) v_temp = _mm512_setzero_si512();

        i++;
    }

    // 4. 处理剩余位 (Remaining bits / Padding)
    // 如果输入耗尽但输出还没填满 (i < rlen)
    if (i < rlen)
    {
        // 计算最后剩余的部分: r[i] = (temp >> shr_pos) & mask
        __m128i k_shr = _mm_cvtsi32_si128(shr_pos);
        __m512i v_final = _mm512_srl_epi64(v_temp, k_shr);
        
        r[i] = _mm512_and_si512(v_final, v_mask);
        i++;
    }

    // 5. 填充高位零
    for (; i < rlen; i++) {
        r[i] = _mm512_setzero_si512();
    }
}

void avx1_withIFMA(pip_ctx *ctx, int eff_times)
{
	int i;
	buf1_st *dblp;
	ep_t tmp;
    ep_t *p1, *p2;
    point52_t avx1_out[AVX_WAY];
    uint64_t u1x[AVX_WAY][HT_NWORDS], u1y[AVX_WAY][HT_NWORDS], u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    
    struct timespec t1, t2;

    /*
    // --- 阶段1: 数据转换 --- 
    clock_gettime(CLOCK_MONOTONIC, &t1);
	// once simdPADDaff
    
	for(i = 0; i < eff_times; i++)
	{
		dblp = &ctx->buf1[i];
        p1 = (ep_t *)(dblp->P0);
        p2 = (ep_t *)(dblp->P1);
        mpi_conv_64to52(u1x[i], (*p1)->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u1y[i], (*p1)->y, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2x[i], (*p2)->x, HT_NWORDS, 6); // convert to radix-52
        mpi_conv_64to52(u2y[i], (*p2)->y, HT_NWORDS, 6); // convert to radix-52
		// ep_add_projc(tmp, *(dblp->P0), *(dblp->P1));
	}
    clock_gettime(CLOCK_MONOTONIC, &t2);
    total_avx_prepare += get_elapsed_time(t1, t2);
    */
   //--- 阶段1: 数据转换 (AVX-512 垂直并行优化版) --- /
    clock_gettime(CLOCK_MONOTONIC, &t1);

    // 1. 定义垂直并行的输入/输出缓冲区
    // htfe_t 是 __m512i[HT_NWORDS]，但在输入时我们只需要用到前 INPUT_LIMBS 个向量
    htfe_t in_u1x, in_u1y, in_u2x, in_u2y;
    htfe_t out_u1x, out_u1y, out_u2x, out_u2y;

    // 临时数组，用于构建垂直向量 (8个点的同一位置肢)
    uint64_t t_1x[8], t_1y[8], t_2x[8], t_2y[8];

    // =========================================================
    // 步骤 A: Gather (收集) - 将 AoS 转置为 SoA
    // 我们需要遍历输入的长度 (通常是 6 个 limb)
    // =========================================================
    for (int j = 0; j < 6; j++) 
    {
        // 收集所有 8 个点的第 j 个 limb
        for (int k = 0; k < eff_times; k++) {
            dblp = &ctx->buf1[k];
            // 注意：这里假设 P0 是指向 ep_t 的指针，ep_t 包含 x, y 数组
            p1 = (ep_t *)(dblp->P0);
            p2 = (ep_t *)(dblp->P1);

            t_1x[k] = (*p1)->x[j];
            t_1y[k] = (*p1)->y[j];
            t_2x[k] = (*p2)->x[j];
            t_2y[k] = (*p2)->y[j];
        }

        // 处理不足 8 个点的情况 (Padding 0)
        // 这一步很重要，避免读取未初始化内存或计算出错误的 NaN
        if (eff_times < 8) {
            int remain = 8 - eff_times;
            memset(&t_1x[eff_times], 0, remain * sizeof(uint64_t));
            memset(&t_1y[eff_times], 0, remain * sizeof(uint64_t));
            memset(&t_2x[eff_times], 0, remain * sizeof(uint64_t));
            memset(&t_2y[eff_times], 0, remain * sizeof(uint64_t));
        }

        // 加载到 AVX 寄存器 (形成垂直条纹)
        in_u1x[j] = _mm512_loadu_si512(t_1x);
        in_u1y[j] = _mm512_loadu_si512(t_1y);
        in_u2x[j] = _mm512_loadu_si512(t_2x);
        in_u2y[j] = _mm512_loadu_si512(t_2y);
    }
    
    // =========================================================
    // 步骤 B: Convert (转换) - 并行执行格式转换
    // =========================================================
    // 输入长度 6，输出长度 8 (HT_NWORDS)
    mpi_conv_64to52_8x1w(out_u1x, in_u1x, HT_NWORDS , 6);
    mpi_conv_64to52_8x1w(out_u1y, in_u1y, HT_NWORDS , 6);
    mpi_conv_64to52_8x1w(out_u2x, in_u2x, HT_NWORDS , 6);
    mpi_conv_64to52_8x1w(out_u2y, in_u2y, HT_NWORDS , 6);

    // =========================================================
    // 步骤 C: Scatter (分发) - 将 SoA 转置回 AoS
    // 必须遍历输出的长度 (8 个 limb)
    // =========================================================
    // 如果你的 simdPADDaff 还没有修改为支持 htfe_t 输入，需要这步
    // 将结果存回 u1x[eff_times][HT_NWORDS] 数组中
    
    for (int j = 0; j < HT_NWORDS ; j++) 
    {
        uint64_t t_out[8];

        // 1. 处理 u1x
        _mm512_storeu_si512(t_out, out_u1x[j]);
        for(int k = 0; k < eff_times; k++) u1x[k][j] = t_out[k];

        // 2. 处理 u1y
        _mm512_storeu_si512(t_out, out_u1y[j]);
        for(int k = 0; k < eff_times; k++) u1y[k][j] = t_out[k];

        // 3. 处理 u2x
        _mm512_storeu_si512(t_out, out_u2x[j]);
        for(int k = 0; k < eff_times; k++) u2x[k][j] = t_out[k];

        // 4. 处理 u2y
        _mm512_storeu_si512(t_out, out_u2y[j]);
        for(int k = 0; k < eff_times; k++) u2y[k][j] = t_out[k];
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t2);
    total_avx_prepare += get_elapsed_time(t1, t2);
    
    // --- 阶段2: SIMD 核心计算 --- 
    clock_gettime(CLOCK_MONOTONIC, &t1);

    simdPADDaff(u1x, u1y, u2x, u2y, avx1_out);

    clock_gettime(CLOCK_MONOTONIC, &t2);
    total_avx_calc += get_elapsed_time(t1, t2);	

    /* --- 阶段3: 结果入队 --- */
    clock_gettime(CLOCK_MONOTONIC, &t1);
    for(i = 0; i < eff_times; i++)
    {
        dblp = &ctx->buf1[i];
        if(point52_is_infty(ctx->bucket[dblp->id].B))
        {
            point52_copy(ctx->bucket[dblp->id].B, avx1_out[i]);
        }
        else
        {
            queue_in(&(ctx->T), avx1_out[i], dblp->id);
        }
    }

	ctx->cnt1 = 0;
    clock_gettime(CLOCK_MONOTONIC, &t2);
    total_avx_queue += get_elapsed_time(t1, t2);

    /* --- 阶段4: 状态转移 (重点测试) --- */
    clock_gettime(CLOCK_MONOTONIC, &t1);

	/* deal with the avx1_out, change state in T*/
	state_trans_T_withIFMA(ctx);
    
    clock_gettime(CLOCK_MONOTONIC, &t2);
    total_avx_trans += get_elapsed_time(t1, t2);
}


double total_time_empty=0.0;
double total_time_fill=0.0;
double total_time_avx=0.0;

void put_into_bucket_withIFMA(const ep_t *P, const uint32_t b_id, int p_id, pip_ctx *ctx)
{
	/* judge if the bucket wait is empty*/
    struct timespec t_start, t_end;
    
    // 开始计时
    clock_gettime(CLOCK_MONOTONIC, &t_start);
	if(ctx->bucket[b_id].wait == NULL)
	{

		// printf("%d ", p_id);
		ctx->bucket[b_id].wait = &(P[p_id]);
        
        clock_gettime(CLOCK_MONOTONIC, &t_end);
        total_time_empty += get_elapsed_time(t_start, t_end);
	}

	else
	{
        struct timespec t_mid;

		ctx->buf1[ctx->cnt1].id = b_id;
		ctx->buf1[ctx->cnt1].P0 = ctx->bucket[b_id].wait;
		ctx->buf1[ctx->cnt1].P1 = &(P[p_id]);
		ctx->cnt1++;
		ctx->bucket[b_id].wait = NULL;

        clock_gettime(CLOCK_MONOTONIC, &t_mid); // 逻辑处理结束，AVX开始前
        total_time_fill += get_elapsed_time(t_start, t_mid);
		/* judge if buf1 is full*/
		if(ctx->cnt1 == AVX_WAY)
		{
            struct timespec avx_start, avx_end;
            clock_gettime(CLOCK_MONOTONIC, &avx_start);
			/* avx1_in is full, then do avx1*/
			avx1_withIFMA(ctx, AVX_WAY);

            clock_gettime(CLOCK_MONOTONIC, &avx_end);
            total_time_avx += get_elapsed_time(avx_start, avx_end);
		}
	}
}
void mont64_to_mont52(point52_t q, ep_t p)
{
    ep_t tmp;
    fp_t mont_a;
    mont_a[0] = 0x44f6480ea8e9b9af;
    mont_a[1] = 0xa96f7d65766c8fe4;
    mont_a[2] = 0xe82efd4228b540fe;
    mont_a[3] = 0x6723e5f0ade53b2e;
    mont_a[4] = 0x25ff6eb6fdd4230a;
    mont_a[5] = 0x14c8ee06ef23c24a;

    ep_copy(tmp, p);
    fp_mul(tmp->x, tmp->x, mont_a);
    fp_mul(tmp->y, tmp->y, mont_a);
    fp_mul(tmp->z, tmp->z, mont_a);
    mpi_conv_64to52(q->x, tmp->x, HT_NWORDS, 6);
    mpi_conv_64to52(q->y, tmp->y, HT_NWORDS, 6);
    mpi_conv_64to52(q->z, tmp->z, HT_NWORDS, 6);
}

void mont52_to_mont64(point52_t q, ep_t p)
{
    ep_t tmp;
    fp_t mont_b;

    mont_b[0] = 0x0;
    mont_b[1] = 0x0;
    mont_b[2] = 0x0;
    mont_b[3] = 0x0;
    mont_b[4] = 0x0;
    mont_b[5] = 0x100000000;
    ep_set_infty(tmp);

    mpi_conv_52to64(tmp->x, q->x, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->y, q->y, 6, HT_NWORDS);
    mpi_conv_52to64(tmp->z, q->z, 6, HT_NWORDS);
    fp_mul(p->x, tmp->x, mont_b);
    fp_mul(p->y, tmp->y, mont_b);
    fp_mul(p->z, tmp->z, mont_b); 
    // if(!ep_is_infty(p))
    // {
    //     p->coord = 2;
    // }
    // else
    // {
    //     p->coord = BASIC;
    // }

    // if(!ep_on_curve(p))
    // {
    //     printf("on curve error\n");
    // }
}



void pippenger_ifma(const bn_t *k, const ep_t *P, int num, ep_t result, ep_t *M, pip_ctx *ctx)
{
    int WBITS = ctx->WBITS;
    int NWINS = ctx->NWINS;
    int BUCKETNUM = ctx->BUCKETNUM;

    ep_t Gu[NWINS], Tu[NWINS];
    uint32_t b;
    int i, a;
    uint64_t u2x[AVX_WAY][HT_NWORDS], u2y[AVX_WAY][HT_NWORDS];
    ep_t tmp;

    double t_flush_avx1 = 0.0;
    double t_flush_avx2 = 0.0;
    double t_flush_buckets = 0.0;
    struct timespec ts_start, ts_end; 

    double time_accumulation = 0.0;

    for (i = 0; i < NWINS; i++)
    {
        // printf("第%d列 \n", i);
        init(ctx);
        for (int j = 0; j < num; j++)
        {
            b = get_wval(k[j], i, ctx->WBITS);
            if (b != 0)
            {
                // ep_add_projc(ctx->bucket[(uint32_t)b - 1].B, ctx->bucket[(uint32_t)b - 1].B, P[j]);
                put_into_bucket_withIFMA(P, (uint32_t)b - 1, j, ctx);
            }
        }

        /* deal with the remains in buf*/
		/* remain buf1*/
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

		if(ctx->cnt1 != 0)
		{
			avx1_withIFMA(ctx, ctx->cnt1);
		}
        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        t_flush_avx1 += get_elapsed_time(ts_start, ts_end);

		/* remain buf2*/
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

		while((!queue_is_empty(&(ctx->T))) || (ctx->cnt2 != 0))
		{
			avx2_withIFMA(ctx, ctx->cnt2);
			state_trans_T_withIFMA(ctx);
		}

        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        t_flush_avx2 += get_elapsed_time(ts_start, ts_end);

        /* deal with the remains in buf0*/
        clock_gettime(CLOCK_MONOTONIC, &ts_start);

        for (a = 0; a < BUCKETNUM; a++)
        {
            if (ctx->bucket[a].wait != NULL)
            {
                // ep_add_projc(ctx->bucket[ctx->dblpoint[a].id].B, ctx->bucket[ctx->dblpoint[a].id].B, *(ep_t *)(ctx->dblpoint[a].P0));
                if(point52_is_infty(ctx->bucket[a].B))
                {
                    mont64_to_mont52(ctx->bucket[a].B, *(ctx->bucket[a].wait));
                }
                else
                {
                    ctx->avx3_out[ctx->avx3_cnt] = &(ctx->bucket[a].B);
                    mpi_conv_64to52(u2x[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->x, HT_NWORDS, 6); // convert to radix-52
                    mpi_conv_64to52(u2y[ctx->avx3_cnt], (*(ctx->bucket[a].wait))->y, HT_NWORDS, 6); // convert to radix-52
                    // mpi_conv_64to52(u2z[ctx->avx3_cnt], (*(ep_t *)(ctx->dblpoint[a].P0))->z, HT_NWORDS, 6); // convert to radix-52
                    ctx->avx3_cnt++;

                    if(ctx->avx3_cnt == AVX_WAY)
                    {
                        simdPADDmix(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
                        ctx->avx3_cnt = 0;
                    }
                }
            }
        }
        if(ctx->avx3_cnt)
        {
            simdPADDmix(ctx->avx3_out, u2x, u2y, ctx->avx3_cnt);
        }

        clock_gettime(CLOCK_MONOTONIC, &ts_end);
        t_flush_buckets += get_elapsed_time(ts_start, ts_end);


        /* Summary points in buckets*/
        struct timespec t_accum_start, t_accum_end;
        clock_gettime(CLOCK_MONOTONIC, &t_accum_start);

        mont52_to_mont64(ctx->bucket[BUCKETNUM - 1].B, M[0]);
        M[0]->coord = PROJC;
        ep_copy(Gu[i], M[0]);
        tmp->coord = PROJC;

        for (int j = 1; j < BUCKETNUM; j++)
        {
            mont52_to_mont64(ctx->bucket[BUCKETNUM - 1 - j].B, tmp);
            ep_add_projc(M[j], tmp, M[j - 1]);
            ep_add_projc(Gu[i], Gu[i], M[j]);
        }

        clock_gettime(CLOCK_MONOTONIC, &t_accum_end); // --- 计时结束 ---
        
        // 2. 计算并打印时间
        time_accumulation += get_elapsed_time(t_accum_start, t_accum_end);
        
        // printf("第%d列 \n", i);
        // ep_norm(Gu[i], Gu[i]);
        // ep_print(Gu[i]);
        // if(Gu[i]->coord != PROJC)
        // {
        //     printf("Gu[i]->coord = \n", Gu[i]->coord);
        // }
        // if(!ep_on_curve(Gu[i]))
        // {
        //     printf("% d on curve error\n", i);
        // }



    }


    printf("------- 详细性能拆解 -------\n");
    printf(">> put_into_bucket 内部:\n");
    printf("   - Bucket为空直接赋值 : %.6f s\n", total_time_empty);
    printf("   - 填充Buffer逻辑     : %.6f s\n", total_time_fill);
    printf("   - AVX函数总调用耗时  : %.6f s\n", total_time_avx);
    
    printf(">> avx1_withIFMA 内部细分:\n");
    printf("   - [1] 数据转换(Prepare): %.6f s\n", total_avx_prepare);
    printf("   - [2] SIMD计算(Calc)   : %.6f s\n", total_avx_calc);
    printf("   - [3] 结果入队(Queue)  : %.6f s\n", total_avx_queue);
    printf("   - [4] 状态转移(Trans)  : %.6f s\n", total_avx_trans);
    printf("---------------------------\n");

    printf(">> 状态转移内部 (state_trans):\n");
    printf("   - AVX2计算累计耗时 : %.6f s\n", total_time_avx2);
    printf("---------------------------\n");

    printf(">> 收尾阶段耗时 (Cleanup Phase):\n");
    printf("   - Flush buf1 (AVX1)    : %.6f s\n", t_flush_avx1);
    printf("   - Flush buf2 (AVX2/Q)  : %.6f s\n", t_flush_avx2);
    printf("   - Flush Buckets (Loop) : %.6f s\n", t_flush_buckets);
    printf("---------------------------\n");

    printf(">> 桶累加耗时 (Accumulation): %.6f s\n", time_accumulation);


    struct timespec t_final_start, t_final_end;
    double time_final_combine = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &t_final_start);

    ep_copy(Tu[0], Gu[NWINS - 1]);
    for (int i = 1; i < NWINS; i++)
    {
        for (int j = 0; j < WBITS; j++)
        {
            ep_dbl_projc(Tu[i - 1], Tu[i - 1]);
        }
        ep_add_projc(Tu[i], Gu[NWINS - 1 - i], Tu[i - 1]);
    }
    ep_copy(result, Tu[NWINS - 1]);

    clock_gettime(CLOCK_MONOTONIC, &t_final_end); 

    time_final_combine = get_elapsed_time(t_final_start, t_final_end);
    printf(">> 最终窗口合并耗时 : %.6f 秒\n", time_final_combine);
}