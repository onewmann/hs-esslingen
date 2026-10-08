function test_phase_correlation()
%TEST_PHASE_CORRELATION  Reference shift between two frames.
    seed_rng(5);
    [~, fr] = simulate_pulsed_image(128, [3 -2], 1);
    [dx, dy] = phase_correlation(fr(:, :, 1), fr(:, :, 2));
    assert(abs(dx - 3) < 0.05 && abs(dy + 2) < 0.05, sprintf('integer shift: %.3f %.3f', dx, dy));
    [~, fr] = simulate_pulsed_image(128, [-1.5 0.75], 1);
    [dx, dy] = phase_correlation(fr(:, :, 1), fr(:, :, 3));
    assert(abs(dx + 3) < 0.15 && abs(dy - 1.5) < 0.15, sprintf('sub-pixel shift: %.3f %.3f', dx, dy));
end
