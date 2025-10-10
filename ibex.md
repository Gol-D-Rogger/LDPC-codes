Step 1: Set likelihood_ram_data_in base two soft bits and calculate the initial syndrome cn base one hard bit. If syndrome are all zero, then decode pass; otherwise, go to step 2.
Step 2: Calculate syndrome_weight and set likelihool_level.
Setp 3: Iteration

for (iter = 0; iter < max_iter; iter++) 
    for (i = 0; i < N; i++) 
    Step 3.1: Update likelihood_new[i]: The likelihood of the current variable node hard decision value is error
        casex ({likelihood_init, likelihood_old[1:0]}): In the first iteration, likelihood_init=1, initialize likelihood of the current variable node.
            3'b0xx: likelihood = likelihood_old;
            3'b100: likelihood = likelihood_level[0];
            3'b101: likelihood = likelihood_level[1];
            3'b110: likelihood = likelihood_level[2];
            3'b111: likelihood = likelihood_level[3];
            default: likelihood = likelihood_old;
        endcase
        a = sum(syndrome[i]): The sum of the syndrome connected to the current variable node.
        clamp = (likelihood == likelihood_level[0])
        flapped_in = likelihood_old >= LIKELIHOOD_FLIP_THR;
        casex({be_aggressive, flipped_in, a})
            5'bx0000: delta = clamp ? 0 : -1;
            5'bx0001: delta = 0;
            5'b00010: delta = 1;
            5'b00011: delta = 2;
            5'b00100: delta = 3;
            5'b10010: delta = 2;
            5'b10011: delta = 4;
            5'b10100: delta = 6;
            5'bx1000: delta = 0;
            5'bx1001: delta = -1;
            5'bx1010: delta = -2;
            5'bx1011: delta = -3;
            5'bx1100: delta = -4;
            default: delta = 0;
        endcase
        In strobe 1/3/5/7, likelihood_level[0] is the smallest. And the smaller the likelihood, the more reliable the hard decision bits are. So, clamp=1 and sum of syndrome is equal to 0, likelihood_new=likelihood_old. When flipped_in=1, delta is equal to sum of syndrome and tend to filp current hard decision.
        If itr <= last_normal_iter
            likelihood_new[i] = likelihood + delta; // update likelihood
        else // post process
            likelihood_new[i] = likelihood + delta; 
            if prng_post_process[0] && (likelihood_new == LIKELIHOOD_FLIP_THR)
                likelihood_new[i] = LIKELIHOOD_FLIP_THR + 1; // increase the difficulty to flip
            else if prng_post_process[1] && (likelihood_new < LIKELIHOOD_FLIP_THR)
                likelihood_new[i] = LIKELIHOOD_FLIP_THR - 1; // decrease the difficulty to flip
            else if prng_post_process[1] && (likelihood_new == LIKELIHOOD_FLIP_THR - 1) && !flipped_in && a==1
                likelihood_new[i] = LIKELIHOOD_FLIP_THR;(Decrease the difficulty to flip to accelerate the convergence)
            end if
        end if
        flipped_out = likelihood_en & (likelihood_new >= LIKELIHOOD_FLIP_THR);
        toggle = likelihood & (flipped_in ^ flipped_out);

        Step 3.2 update cn
            for g = 1:1:M
                mask_for_row[g] = mask_row[g] ? mask : fade_row[g] ? ~mask : row_enable[g] ? {P{1'b1}} : {P{1'b0}};
                s[g] = {cn[g][P-1-LDPC_MATRIX_DELTA{g}:0], cn[g][P-1:P-LDPC_MATRIX_DELTA{g}]};
                cn[g] = (toggle & mask_for_row[g]) ^ ((state==ST_IDLE))?0:s[g];
            end for

        Step 3.3 Determine if decoding is complete
            if cn == 0
                break;
        end if
    end for
end for

Step 4: Output
flipped_out = [flipped_hi, flipped_lo];
flip_ram_data_out = flip_ram_out_select ? flipped_hi : flipped_lo;
error_count = error_count + sum(flipped_out);
o_data=(error_count==0)?hard_data_r :(hard_data_r ^ flip_ram_data_out);

likelihood_old: likelihood of vn is error of previous iteration, output of likelihood_ram.
likelihood_new: likelihood of vn is error of current iteration.
flipped_in: flip flag of previous iteration.
flipped_out: flip flag of current iteration.
mask_for_row: the mask for circulant.
hard_data_r: raw hard data.
