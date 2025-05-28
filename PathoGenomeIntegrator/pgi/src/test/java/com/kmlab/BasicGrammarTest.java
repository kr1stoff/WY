package com.kmlab;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.Test;

public class BasicGrammarTest {
    @Test
    public void sortTest() {
        int[] totalScores = {3,2,3,4,5,5};
        List<Integer> totalScoreList = new ArrayList<>();
        for (int score : totalScores) {
            totalScoreList.add(score);
        }
        Collections.sort(totalScoreList);
        System.out.println(totalScoreList);
    }
}
