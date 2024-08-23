function [system] = speer()
    eqs = cell(4,1);
    eqs{1} = [-3.94775390625+1.12591552734375i 2 1 1 1; 3.94775390625-1.12591552734375i 2 1 0 0; 3.94775390625-1.12591552734375i 2 0 0 1; 1.31591796875-0.37530517578125i 1 2 1 1; -1.31591796875+0.37530517578125i 1 2 0 0; 1.31591796875-0.37530517578125i 1 1 2 1; 1.31591796875-0.37530517578125i 1 1 1 2; 0.562856865115464-0.365647947415709i 1 1 1 1; 2.6318359375-0.7506103515625i 1 1 1 0; -2.6318359375+0.7506103515625i 1 1 0 1; -0.562856865115464+0.365647947415709i 1 1 0 0; 2.6318359375-0.7506103515625i 1 0 1 1; -1.31591796875+0.37530517578125i 1 0 0 2; -0.562856865115464+0.365647947415709i 1 0 0 1; -1.31591796875+0.37530517578125i 0 2 1 0; -1.31591796875+0.37530517578125i 0 1 2 0; -3.6318359375+0.7506103515625i 0 1 1 1; -0.562856865115464+0.365647947415709i 0 1 1 0; 1 0 1 0 0; -1.31591796875+0.37530517578125i 0 0 2 1; -1.31591796875+0.37530517578125i 0 0 1 2; -0.562856865115464+0.365647947415709i 0 0 1 1; 1 0 0 0 1];
    eqs{2} = [1.31591796875-0.37530517578125i 2 1 1 1; -1.31591796875+0.37530517578125i 2 1 0 0; -1.31591796875+0.37530517578125i 2 0 0 1; 1.31591796875-0.37530517578125i 1 2 1 1; -1.31591796875+0.37530517578125i 1 2 0 0; -3.94775390625+1.12591552734375i 1 1 2 1; 1.31591796875-0.37530517578125i 1 1 1 2; 0.562856865115464-0.365647947415709i 1 1 1 1; 2.6318359375-0.7506103515625i 1 1 1 0; -3.6318359375+0.7506103515625i 1 1 0 1; -0.562856865115464+0.365647947415709i 1 1 0 0; 2.6318359375-0.7506103515625i 1 0 1 1; -1.31591796875+0.37530517578125i 1 0 0 2; -0.562856865115464+0.365647947415709i 1 0 0 1; -1.31591796875+0.37530517578125i 0 2 1 0; 3.94775390625-1.12591552734375i 0 1 2 0; -2.6318359375+0.7506103515625i 0 1 1 1; -0.562856865115464+0.365647947415709i 0 1 1 0; 1 0 1 0 0; 3.94775390625-1.12591552734375i 0 0 2 1; -1.31591796875+0.37530517578125i 0 0 1 2; -0.562856865115464+0.365647947415709i 0 0 1 1; 1 0 0 0 1];
    eqs{3} = [1.31591796875-0.37530517578125i 2 1 1 1; -1.31591796875+0.37530517578125i 2 1 0 0; -1.31591796875+0.37530517578125i 2 0 0 1; -3.94775390625+1.12591552734375i 1 2 1 1; 3.94775390625-1.12591552734375i 1 2 0 0; 1.31591796875-0.37530517578125i 1 1 2 1; 1.31591796875-0.37530517578125i 1 1 1 2; 0.562856865115464-0.365647947415709i 1 1 1 1; -2.6318359375+0.7506103515625i 1 1 1 0; 2.6318359375-0.7506103515625i 1 1 0 1; -0.562856865115464+0.365647947415709i 1 1 0 0; -3.6318359375+0.7506103515625i 1 0 1 1; -1.31591796875+0.37530517578125i 1 0 0 2; -0.562856865115464+0.365647947415709i 1 0 0 1; 1 1 0 0 0; 3.94775390625-1.12591552734375i 0 2 1 0; -1.31591796875+0.37530517578125i 0 1 2 0; 2.6318359375-0.7506103515625i 0 1 1 1; -0.562856865115464+0.365647947415709i 0 1 1 0; -1.31591796875+0.37530517578125i 0 0 2 1; -1.31591796875+0.37530517578125i 0 0 1 2; -0.562856865115464+0.365647947415709i 0 0 1 1; 1 0 0 1 0];
    eqs{4} = [1.31591796875-0.37530517578125i 2 1 1 1; -1.31591796875+0.37530517578125i 2 1 0 0; -1.31591796875+0.37530517578125i 2 0 0 1; 1.31591796875-0.37530517578125i 1 2 1 1; -1.31591796875+0.37530517578125i 1 2 0 0; 1.31591796875-0.37530517578125i 1 1 2 1; -3.94775390625+1.12591552734375i 1 1 1 2; 0.562856865115464-0.365647947415709i 1 1 1 1; -3.6318359375+0.7506103515625i 1 1 1 0; 2.6318359375-0.7506103515625i 1 1 0 1; -0.562856865115464+0.365647947415709i 1 1 0 0; -2.6318359375+0.7506103515625i 1 0 1 1; 3.94775390625-1.12591552734375i 1 0 0 2; -0.562856865115464+0.365647947415709i 1 0 0 1; 1 1 0 0 0; -1.31591796875+0.37530517578125i 0 2 1 0; -1.31591796875+0.37530517578125i 0 1 2 0; 2.6318359375-0.7506103515625i 0 1 1 1; -0.562856865115464+0.365647947415709i 0 1 1 0; -1.31591796875+0.37530517578125i 0 0 2 1; 3.94775390625-1.12591552734375i 0 0 1 2; -0.562856865115464+0.365647947415709i 0 0 1 1; 1 0 0 1 0];

    system = systemstruct(eqs);
end