package com.kmlab.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.junit.Test;

public class RankingCalculatorTest {
    @Test
    public void intMapRankTest() {
        int[] totalScores = {3,2,3,4,5,5};
        List<Integer> totalScoreList = new ArrayList<Integer>();

        for (int score : totalScores) {
            totalScoreList.add(score);
        }

        Map<Integer, Integer> totalScoresRankMap = RankingCalculator.intMapRank(totalScoreList);

        for (Integer key : totalScoresRankMap.keySet()) {
            System.out.println(key + ": " + totalScoresRankMap.get(key));
        }
    }
}
