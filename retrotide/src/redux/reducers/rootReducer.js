import {
  UPDATE_MOLECULE,
  UPDATE_MODULES,
} from '../actions/actions';

const initialState = {
  sessionID: "",
  jobs: [],
  molecule: "",
};

function rootReducer(state = initialState, action) {
  console.log("REDUCER HIT: Action ->", action);

  switch (action.type) {
    case UPDATE_MOLECULE:
      return { ...state, Molecule: action.payload.molecule };

    case UPDATE_MODULES:
      return { ...state, ModuleArray: action.payload.ModuleArray }

    // case SUBMIT_DESIGN:
    //   return { ...state, }
    default:
      return state;
  }
}

export default rootReducer;